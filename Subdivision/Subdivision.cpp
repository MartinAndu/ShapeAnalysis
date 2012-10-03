#include <fstream>
#include <vector>

#include "OBJFileReader.h"
#include "Solid.h"
#include "iterators.h"
#include "Edge.h"
#include "Triangle.h"
#include "SolidDelegate.h"

#define PI 3.14159265

using namespace MeshLib;

double getAlpha(int n){
	double answer;
	if(n>3){
		double center = (0.375 + (0.25 * cos((2.0 * 3.14159265358979) / (double) n)));
		answer = (0.625 - (center * center)) / (double) n;
	} else {
		answer = 3.0 / 16.0;
	}
	return answer;
}

void main(int argc, char *argv[])
{
	// Read in the obj file
	Solid mesh, newMesh;
	newMesh = Solid();
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&mesh, in);
	
	Solid  meshes;

	for (int k=0; k<atoi(argv[3]); k++) {
	
	//if (i==1) newMesh = Solid();

		int numVertices = mesh.numVertices();
		int numFaces = mesh.numFaces();
		int numEdges = mesh.numEdges();
		SolidDelegate delegate;

		// Future new faces
		int *newverts = new int[numFaces*10];	//ABCA'B'C'
		int faceID = 0;

		// Iterate thru all faces
		SolidFaceIterator faces(&mesh);
		for (int newVertsId = numVertices, global = 0; !faces.end(); ++faces,global++) {
			// ---------------------------
			// Get old face vertices
			// ---------------------------
			Solid::tFace face = *faces;
			Solid::tVertex va, vb, vc;
			va = new Vertex();
			vb = new Vertex();
			vc = new Vertex();
			FaceVertexIterator vertices(face);
			for (int i=0; !vertices.end(); ++vertices, i++) {
				Solid::tVertex vertex = *vertices;
				//std::cout << "\r\ng=" << global << ", i = " << i << ", vertex = " << vertex->id();
				if (i==0) va = vertex;
				else if (i==1) vb = vertex;
				else if (i==2) vc = vertex;
			}

			// ---------------------------
			// Calculate new additional vertices
			// ---------------------------
			Vertex *vap, *vbp, *vcp;
				
			vap = delegate.createVertex(&mesh, ++newVertsId);
			vap->point() = (vb->point()+vc->point())/2;

			vbp = delegate.createVertex(&mesh, ++newVertsId);
			vbp->point() = (va->point()+vc->point())/2;

			vcp = delegate.createVertex(&mesh, ++newVertsId);
			vcp->point() = (va->point()+vb->point())/2;
		
			// ---------------------------
			// Store new faces' vertices
			// ---------------------------
			newverts[global+(0*numFaces)] = va ->id();
			newverts[global+(1*numFaces)] = vb ->id();
			newverts[global+(2*numFaces)] = vc ->id();
			newverts[global+(3*numFaces)] = vap->id();
			newverts[global+(4*numFaces)] = vbp->id();
			newverts[global+(5*numFaces)] = vcp->id();		
		}

		// ---------------------------
		// Insert ALL vertices into new mesh
		// ---------------------------
		SolidVertexIterator vs(&mesh);
		for (int i=1; !vs.end(); ++vs) {
			Solid::tVertex v = *vs;
			Vertex * vv = delegate.createVertex(&newMesh,  v->id());
			vv->point() = v->point();
			vv->id()    = v->id();
	
		}
	
		// ---------------------------
		// Insert new faces into new mesh
		// ---------------------------
		Face *face1, *face2, *face3, *face4;
		Edge *edge;
		for (int i=0; i<numFaces; i++) {
			//std::cout<<"\r\nInserting face #"<<i;
			int va, vb, vc, vap, vbp, vcp;

			va = newverts[i+(0*numFaces)];
			vb = newverts[i+(1*numFaces)];
			vc = newverts[i+(2*numFaces)];
			vap = newverts[i+(3*numFaces)];
			vbp = newverts[i+(4*numFaces)];
			vcp = newverts[i+(5*numFaces)];

			//  New triangles: AC'B', C'BA', B'A'C, C'A'B'
			//
			//                  A
			//	              /   \
			//               /  1  \
			//              C'______B'
			//             / \  4  / \
			//            / 2 \   / 3 \
			//           B______A'_____C
			//
			// (may be in different order because of normals)

			int* fv1 = new int[3];
			int* fv2 = new int[3];
			int* fv3 = new int[3];
			int* fv4 = new int[3];
		
			fv1[0] = va; // the order of these is important because it will 
			fv1[1] = vcp;// define the normal of the triangle.
			fv1[2] = vbp;
			face1 = delegate.createFace(&newMesh, fv1, i+(numFaces*0));
		
			fv2[0] = vcp;
			fv2[1] = vb;
			fv2[2] = vap;
			face2 = delegate.createFace(&newMesh, fv2, i+(numFaces*1));

			fv3[0] = vbp;
			fv3[1] = vap;
			fv3[2] = vc;
			face3 = delegate.createFace(&newMesh, fv3, i+(numFaces*2));

			fv4[0] = vcp;
			fv4[1] = vap;
			fv4[2] = vbp;
			face4 = delegate.createFace(&newMesh, fv4, i+(numFaces*3));
		}
	
	
		// Subdivide: USE ORIGINAL MESH!! then update using ID
		SolidVertexIterator originalV(&mesh);
		// modify only the original vertices
		for(int i=0; !originalV.end() && i<numVertices; i++, ++originalV)  {	
			Solid::tVertex vert = *originalV;

			// 1. sum neighbors
			Point summation = Point();
			int nn=0;	//n neighbors
			MeshLib::VertexVertexIterator neighbors(vert);
			for (Solid::tVertex vn = *neighbors; !neighbors.end(); nn++, ++neighbors) summation+=vn->point();
			double n = (double) nn;
			// 2. calculate constant
			//double s = (1/n) * (5/8 - ( (3/8 + ( ( (double) cos ((2*PI) / n) )/4 ) ) * 2) ) ;
			double beta = (n>3) ? /*(3/(8*n))*/getAlpha(n) : (3/16);
			std::cout<<"\r\n i = "<<i<< " n= "<<n<< " beta= "<<beta<<"\t sum= "<<summation(0)<<", "<<summation(1)<<", "<<summation(2);
			
			// 3. calculate and update with new position
			Point old = vert->point();
			Point newpoint = old*(1-(n*beta)) + summation*beta;	
			vert->point() = newpoint;
			newMesh.idVertex(vert->id())->point() = newpoint;	
		}
	std::cout<<"\r\n ================================ \r\n Finished pass "<<k<<"\r\n============================";
		newMesh.copy(mesh);
	//	mesh = newMesh;
	}

	// Write out the resultant obj file
	int vObjID = 1;
	std::map<int, int> vidToObjID;

	std::ofstream os(argv[2]);
	
	SolidVertexIterator iter(&newMesh);
	
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
	
	SolidFaceIterator fiter(&newMesh);
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
