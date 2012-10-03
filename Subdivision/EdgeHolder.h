#pragma once
class EdgeHolder
{
public:
	int v1, v2, vNew;

	EdgeHolder(int v11, int v22, int vNeww) {v1=v11;v2=v22;vNew=vNeww};
	~EdgeHolder(void);
};

