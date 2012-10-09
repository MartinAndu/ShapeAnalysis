#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "Edge.h"
#include "Triangle.h"
#include "SolidDelegate.h"


#define PI 3.14159265

#define IDX(i,j) (((i<j)?i:j)*nVert)+((i<j)?j:i)

#define VA 0
#define VB 1
#define VC 2
#define VAP 3
#define VBP 4
#define VCP 5



using namespace MeshLib;

double* alphas;

double getAlpha(int n){
	double answer;
	if(n>3){
		double center = (0.375 + (0.25 * cos((2.0 * 3.14159265358979) / (double) n)));
		answer = (0.625 - (center * center)) / (double) n;
	} else {
		answer = 3.0 / 16.0;
		//if (n<3) std::cout<<"\r\nBoundary";
	}
	return answer;
}

SolidDelegate delegate;

// refine mesh
Solid* firstPass(Solid* original) {

	Solid *mesh = new Solid();

	Solid *refinedMesh = new Solid();

	original->copyto(*mesh);

	int nEdge = mesh->numEdges();
	int nVert = mesh->numVertices();
	int nFace = mesh->numFaces();

	int newVertIds =  nVert;
	int newFaceIds = nFace;

	int *newFaceVertIds; // will store 6 vertices IDs PER TRIANGLE.

	newFaceVertIds = new int[6*nFace];

	int * oppositeVertOfEdge = new int[2*(6*nVert)];


	// Lookup table for marked edges.
	// Contains the IDs of the odd (new) vertices for each edge.
	// This could be done with an array of EdgeKeys, but it would have a 
	// linear access time (having to lookup the whole array each time),
	// so we sacrifice some space to get better times.
	int *marked = new int[nVert*nVert];
	for (int i=0; i<nVert*nVert; i++) marked[i] = -1;
	
	int currentFace = 0;
	for (SolidFaceIterator F(original); !F.end(); ++F, currentFace++) {
		Solid::tFace f = *F;
		// 1ST STEP: CREATE NEW VERTICES (SPLIT EDGES).
			
		Solid::tVertex vap, vbp, vcp, va, vb, vc;	// all 6 points for 4 tris

		Solid::tEdge e1, e2, e3;
		Solid::tVertex v31, v32, v33;

		// intermezzo: get the original A,B,C
		int k = 0;
		for (FaceVertexIterator V(f); !V.end(); ++V, k++) {
			Solid::tVertex vv = *V;
			if		(k==0) va = vv;
			else if	(k==1) vb = vv;
			else if (k==2) vc = vv;
		}

		k = 0;
		for(FaceEdgeIterator E(f); !E.end(); ++E, k++) {
			Solid::tEdge e = *E;
			Solid::tVertex v1,v2,v3;
			e->get_vertices(v1,v2);

			// At this point we have an edge and its two
			// vertices, say A and B. We want to define
			// a new vertex V3 between those two.
			//
			//  V1---------==---------V2

			int v3tentative = marked[IDX(v1->id(),v2->id())];

			/*for (int i=0; i<nVert; i++) {
				for (int j=0; j<nVert; j++) 
				if (marked[IDX(i,j)]!=-1) std::cout<<"\r\nPOS("<<i<<","<<j<<") = "<<IDX(i,j)<<" = "<<marked[IDX(i,j)];
			}*/
	


			// we are analyzing edge v1,v2 now, so we want to find vOp

			int vOp = -1;

			if ((va->id() != v1->id()) && (va->id() != v2->id())) vOp = va->id();
			else if ((vb->id() != v1->id()) && (vb->id() != v2->id())) vOp = vb->id();
			else if ((vc->id() != v1->id()) && (vc->id() != v2->id())) vOp = vc->id();

			assert(vOp>0);

			//------------------------------------------------
			// Determine v3
			//------------------------------------------------
			//
			// Each v3 found will be one of vap, vbp, vcp.
			//
			//------------------------------------------------

			// in this case, the edge was already processed by another face.
			if (v3tentative>0) {
				v3 = mesh->idVertex(v3tentative);

				// we will take this opportunity to store one of 
				// the opposite even vertices to this new odd vertex,
				// the other one was already discovered since the 
				// split already happened
				//
				oppositeVertOfEdge[v3->id()+(6*nVert)] = vOp;
			}
			// in this case, the edge was not yet processed. We need to create it, split the edge etc.
			else {
				// first: split edge, create vertex	
				v3 = delegate.edgeSplit(mesh, mesh->idEdge(v1->id(),v2->id()));
				v3->point() = ((v1->point()+v2->point())*(0.252));	// juicy part
				v3->id() = ++newVertIds;			
				// second: update lookup table
				assert(marked[IDX(v1->id(),v2->id())] < 0);
				marked[IDX(v1->id(),v2->id())]= v3->id();

				// opposite of v3, used to update odd vertices
				oppositeVertOfEdge[v3->id()] = vOp;
			}

			// at this point we have successfully splitted an edge or,
			// if it was already split, we have recovered the odd vertex
			// from the mesh. The only thing left is to assign v3, the odd
			// vertex, to one of the following: vap, vbp, vcp. We will use
			// the following convention:
			//
			// vap = 1st, vbp = 2nd, vcp = 3rd iterations, i.e. k=0,1,2.
			if		(k==0) vbp = v3;
			else if	(k==1) vcp = v3;
			else if (k==2) vap = v3;

			

		}		

	


		// Here we have all 6 vertices. Store them.

		newFaceVertIds[currentFace+(VA*nFace)] = va-> id();
		newFaceVertIds[currentFace+(VB*nFace)] = vb-> id();
		newFaceVertIds[currentFace+(VC*nFace)] = vc-> id();
		newFaceVertIds[currentFace+(VAP*nFace)] = vap->id();
		newFaceVertIds[currentFace+(VBP*nFace)] = vbp->id();
		newFaceVertIds[currentFace+(VCP*nFace)] = vcp->id();
	}


	// Intermezzo: copy all vertices to a fresh mesh
	int j=0;
	for(SolidVertexIterator vi(mesh); !vi.end(); ++vi, j++) {
		Solid::tVertex origV = *vi;
		Vertex *genV = delegate.createVertex(refinedMesh, origV->id());
		genV->point() = origV->point();
		genV->id() = origV->id();
		//genV->boundary() = true;
	}


	//------------------------------------------------
	// 2ND STEP: CREATE NEW EDGES AND FACES.
	//------------------------------------------------
	//
	// At this point all the new (odd) vertices were already created. 
	//
	// We have a mesh whose triangles look something like this:
	//
	//
	//                  A
	//	              /   \
	//               /     \
	//              C'      B'
	//             /         \
	//            /           \
	//           B______A'_____C
	//
	//				Figure 1
	//		An individual triangle with
	//		  even and odd vertices
	//
	// 
	// We want to split each triangle ABC into 4 triangles, namely
	// AC'B', BC'A', A'B'C and A'B'C' (may be in different order 
	// because of normals).
	//
	// Let us examine carefully the mesh we have right now.
	// It will look something like this (projected on a 2D plane, naturally)
	// 
	//						...
	//				
	// 
	//                  E------O-------E          
	//	              /   \          /   \        ...
	//               /     \   2    /     \      
	//              O       O      O       O     
	//             /    1    \    /    3    \  
	//            /           \  /           \ 
	//            E------O------E------O------E             
	//	        /   \         /   \          /   
	//         /     \   6   /     \   4    /   
	//        O       O     O       O      O    
	//       /    7    \   /    5    \    /    ...
	//      /           \ /           \  /    
	//     E------O------E------O-------E
	//
	//						Figure 2
	//		The modified mesh with odd vertices
	//
	// 
	// It is the case here that we may well reconstruct the entire mesh
	// from scratch, already with the new edges/faces between the odd
	// vertices, just with the vertex information we already have at this 
	// point, and the topology will in fact be preserved. 
	//
	// For example, take triangles 1 and 2. We want a total of 8 triangles,
	// 4 for each. They would share 2 even vertices and 1 odd vertex.
	// If we create them independently,	the "central"  triangles (6 and 7 in
	// figure 3) will end up "naturally" sharing the odd vertex, thus 
	// preserving connectivity. 
	//
	//
	//                 E------O------E
	//	              /   \ 5'/ \ 7'/
	//               /  1' \ / 6'\ / 
	//              O-------O-----O  
	//             / \  4' / \ 8'/
	//            / 2'\   / 3'\ /
	//           E------O------E
	//		
	//					Figure 3 
	//	    close up of triangles 6 and 7 
	//       from previous figure after
	//        going through refinement
	//
	// It is crucial to understand that we have no guarantees as to the
	// specific order we will be creating the new triangles i' , i from 
	// 0 to (4*nFace)-1. Because of that it is important that our method
	// is robust enough to preserve correct topology no matter what order
	// we create the new triangles.
	//
	//------------------------------------------------

	// iterate thru all triangles, each with 6 vertices.
	int splitFaceId = 0;
	for(int i=0; i<nFace; i++) {

		int va, vb, vc, vap, vbp, vcp;
		Solid::tVertex vva, vvb, vvc, vvap, vvbp, vvcp, op1_vap, op2_vap, op1_vbp, op2_vbp, op1_vcp, op2_vcp;

		va  = newFaceVertIds[i+VA *nFace];
		vb  = newFaceVertIds[i+VB *nFace];
		vc  = newFaceVertIds[i+VC *nFace];
		vap = newFaceVertIds[i+VAP*nFace];
		vbp = newFaceVertIds[i+VBP*nFace];
		vcp = newFaceVertIds[i+VCP*nFace];

		int tri_A_Cp_Bp [3] = { va, vcp, vbp }; 
		int tri_Cp_B_Ap [3] = { vcp, vb, vap }; 
		int tri_Bp_Ap_C [3] = { vbp, vap, vc }; 
		int tri_Cp_Ap_Bp[3] = { vcp, vap, vbp }; 



		//                  A
		//	              /   \
		//               /  1  \
		//              C'______B'
		//             / \  4  / \
		//            / 2 \   / 3 \
		//           B______A'_____C
		
		// we can now update the vertices with the new coordinates (phase 2)
		//
		// ODDS: VAp, VBp, VCp - UPDATE FIRST!
		// EVEN: VA, VB, VC - UPDATE SECOND!
		//
		vva  = refinedMesh->idVertex(va);
		vvb  = refinedMesh->idVertex(vb);
		vvc  = refinedMesh->idVertex(vc);
		vvap = refinedMesh->idVertex(vap);
		vvbp = refinedMesh->idVertex(vbp);
		vvcp = refinedMesh->idVertex(vcp);

		op1_vap = refinedMesh->idVertex(oppositeVertOfEdge[vap]);
		op2_vap = refinedMesh->idVertex(oppositeVertOfEdge[vap+ 6*nVert]);
		assert(op1_vap!=op2_vap); // do not accept degenerate cases
		op1_vbp = refinedMesh->idVertex(oppositeVertOfEdge[vbp]);
		op2_vbp = refinedMesh->idVertex(oppositeVertOfEdge[vbp + 6*nVert]);
		assert(op1_vbp!=op2_vbp); // do not accept degenerate cases
		op1_vcp = refinedMesh->idVertex(oppositeVertOfEdge[vcp]);
		op2_vcp = refinedMesh->idVertex(oppositeVertOfEdge[vcp + 6*nVert]);
		assert(op1_vcp!=op2_vcp); // do not accept degenerate cases
		

		// UPDATING ODDS
		//
		//        A = 1/8
		//       / \
		//      /   \
		//     /     \
		//    /       \
		// = B----V----C = 3/8
		//3/8 \       /
		//     \     /
		//      \   /
		//       \ /
		//        D = 1/8
		//
		// V = 0.125*(A+D) + 0.375*(C+B)
		//
		// Currently we have V = 0.375*(C+B). We need to add
		// V+=0.125*(A+D) to get the correct point
		//0.0635
		// 
	
		//std::cout<<"\r\n\r\n\r\nAP, \t"<<op1_vap->point()[0]<<"\t\t"<<op2_vap->point()[1];
		//std::cout<<"\r\nBP, \t"<<op1_vbp->point()[0]<<"\t\t"<<op2_vbp->point()[1];
		//std::cout<<"\r\nCP, \t"<<op1_vcp->point()[0]<<"\t\t"<<op2_vcp->point()[1];
		vvap->point() += ((op1_vap->point()+op2_vap->point())*0.125);
		vvbp->point() += ((op1_vbp->point()+op2_vbp->point())*0.125);
		vvcp->point() += ((op1_vcp->point()+op2_vcp->point())*0.125);
		
		
		


		// we can now create 4 triangles:
		delegate.createFace(refinedMesh, tri_A_Cp_Bp,  splitFaceId++);
		delegate.createFace(refinedMesh, tri_Cp_B_Ap,  splitFaceId++);
		delegate.createFace(refinedMesh, tri_Bp_Ap_C,  splitFaceId++);
		delegate.createFace(refinedMesh, tri_Cp_Ap_Bp, splitFaceId++);
	}
	
	// all triangles created; update EVEN verts
	// use REFINEDMESH, because the new neighbors will be there.
	SolidVertexIterator originalV(refinedMesh);
	for(int i=0; !originalV.end() && i<nVert; i++, ++originalV)  {	
		Solid::tVertex vert = *originalV;

		//=======================
		// 1. sum neighbors
		//=======================
		Point summation = Point();
		int nn=0;	//n neighbors
		MeshLib::VertexVertexIterator neighbors(vert);
		for (Solid::tVertex vn = *neighbors; !neighbors.end(); nn++, ++neighbors) summation+=vn->point();
		double n = (double) nn;
		
		//=======================
		// 2. calculate constant
		//=======================
		double beta = (n>3) ? alphas[nn] : (3/16);
		
		//std::cout<<"\r\n i = "<<i<< " n= "<<n<< " beta= "<<beta<<"\t sum= "<<summation(0)<<", "<<summation(1)<<", "<<summation(2);
			
		//=======================
		// 3. calculate and update with new position
		//=======================
		Point old = vert->point();
		Point newpoint = old*(1-(beta*n)) + summation*beta;	
		//if (n==2) newpoint = old*0.75 + summation*0.125;	//crease/boundary
		refinedMesh->idVertex(vert->id())->point() = newpoint;
	}
	
	
	return refinedMesh;
}



void main(int argc, char *argv[])
{
	alphas = new double[30];
	for (int i=0;i<30;i++) alphas[i]=getAlpha(i);

	// Read in the obj file
	Solid mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);

	for (int i=0; i<atoi(argv[3]); i++) {
		mesh = *firstPass(&mesh);
		std::cout<<"\r\nFinished pass #"<<i+1<<", working..."; 
	}
	std::cout<<"\r\n\r\nDone!";



	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);
	
	SolidVertexIterator iter(&mesh);
	
	for(int i=0; !iter.end(); ++iter, i++)
	{
		Vertex *v = *iter;
		Point p = v->point();
		os << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		//std::cout << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
		vidToObjID[v->id()] = vObjID++;
	}
	os << "# " << (unsigned int)mesh.numVertices() << " vertices" << std::endl;

	float u = 0.0, v = 0.0;
	for(iter.reset(); !iter.end(); ++iter)
	{
		Vertex *v = *iter;
		std::string key( "uv" );
		std::string s = Trait::getTraitValue (v->string(), key );
		if( s.length() > 0 )
		{
			sscanf( s.c_str (), "%f %f", &u, &v );
		}
		os << "vt " << 0 << " " << 0 << std::endl;
	}
	os << "# " << (unsigned int)mesh.numVertices() << " texture coordinates" << std::endl;
	
	SolidFaceIterator fiter(&mesh);
	for(; !fiter.end(); ++fiter)
	{
		Face *f = *fiter;
		FaceVertexIterator viter(f);
		os << "f " ;
		for(; !viter.end(); ++viter)
		{
			Vertex *v = *viter;
			os << vidToObjID[v->id()] << "/" << vidToObjID[v->id()] << " ";
		}
		os << std::endl;
	}
	os.close();
}
