typedef struct Vec3f {
	float x;
	float y;
	float z;

} Vec3f;

__kernel void MarchingCubeCL(__global float* data, __global Vec3f* xyz , __global Vec3f* pos, float r, int N,int pN) {

	int id = get_global_id(0);
	int i = 0;
	//printf("%d : %d", id, N);

	if (id < N) {
		for (i = 0; i < pN; i++) {
			data[id] += (r * r) / ((xyz[id].x - pos[i].x)*(xyz[id].x - pos[i].x) + (xyz[id].y - pos[i].y)*(xyz[id].y - pos[i].y) + (xyz[id].z - pos[i].z)*(xyz[id].z - pos[i].z));
		}
		//printf("%f", data[id]);
	}
}

__kernel void zeroBuffer(__global float* data, int N) {
	int id = get_global_id(0);

	if(id < N)data[id] = 0;
}

