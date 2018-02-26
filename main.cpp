#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/colon.h>
#include <igl/slice.h>
#include <igl/upsample.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/harmonic.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>

using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;

void printHelpExit() {
	printf("Invalid command line arguments specified!\n\n");

	printf("USAGE: hole_fixer [options]\n\n");

	printf("OPTIONS: \n");

	printf("  -in\t\t\ttarget mesh file in .off-format, with a hole\n");
	printf("  -out\t\t\toutput mesh file in .off-format\n");
	printf("  -outfaces\t\tHow many faces to decimate the mesh to\n");
	printf("  -upsample\t\tHow much upsampling to use when creating the patch\n");

	exit(1);
}

const char* findToken(const char* param, int argc, char* argv[]) {
	const char* token = nullptr;
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], param) == 0) {
			if (i + 1 < argc) {
				token = argv[i + 1];
				break;
			}
		}
	}

	if (token == nullptr) {
		printf("Could not find command-line parameter %s\n", param);
		return nullptr;
	}

	return token;
}

const char* parseStringParam(const char* param, int argc, char* argv[]) {
	const char* token = findToken(param, argc, argv);
	return token;
}

bool parseIntParam(const char* param, int argc, char* argv[], unsigned int& out) {
	const char* token = findToken(param, argc, argv);
	if (token == nullptr)
		return false;

	int r = sscanf(token, "%u,", &out);
	if (r != 1 || r == EOF) {
		return false;
	}
	else {
		return true;
	}
}

int main(int argc, char *argv[])
{
	
	//
	// command line parsing.
	//
	const char* inFile = parseStringParam("-in", argc, argv);
	if (inFile == nullptr) printHelpExit();

	const char* outFile = parseStringParam("-out", argc, argv);
	if (outFile == nullptr) printHelpExit();

	unsigned int outFacesN;
	if (!parseIntParam("-outfaces", argc, argv, outFacesN)) printHelpExit();

	unsigned int upsampleN;
	if (!parseIntParam("-upsample", argc, argv, upsampleN)) printHelpExit();
	
	// original mesh vertices and indices. This is the original mesh, which has a hole.
	MatrixXd originalV;
	MatrixXi originalF;

	if (!igl::readOFF(inFile, originalV, originalF)) {
		printHelpExit();
	}

	VectorXi originalLoop; // indices of the boundary of the hole. 
	igl::boundary_loop(originalF, originalLoop);

	if (originalLoop.size() == 0) {
		printf("Mesh has no hole!");
		printHelpExit();
	}

	// upsample the original mesh. this makes fusing the original mesh with the patch much easier.
	igl::upsample(Eigen::MatrixXd(originalV), Eigen::MatrixXi(originalF), originalV, originalF, upsampleN);

	// compute boundary center.
	VectorXd bcenter(3);
	{
		VectorXi R = originalLoop;
		VectorXi C(3); C << 0, 1, 2;
		MatrixXd B;
		MatrixXd A = originalV;
		igl::slice(A, R, C, B);
		bcenter = (1.0f / originalLoop.size()) * B.colwise().sum();
	}

	// a flat patch that fills the hole.
	MatrixXd patchV = MatrixXd(originalLoop.size() + 1, 3); // patch will have an extra vertex for the center vertex.
	MatrixXi patchF = MatrixXi(originalLoop.size(), 3);

	{
		VectorXi R = originalLoop;
		VectorXi C(3); C << 0, 1, 2;
		MatrixXd A = originalV;
		MatrixXd temp1;
		igl::slice(A, R, C, temp1);

		MatrixXd temp2(1, 3);
		temp2 << bcenter(0), bcenter(1), bcenter(2);

		// patch vertices will be the boundary vertices, plus a center vertex. concat these together.
		igl::cat(1, temp1, temp2, patchV);

		// create triangles that connect boundary vertices and the center vertex.
		for (int i = 0; i < originalLoop.size(); ++i) {
			patchF(i, 2) = (int)originalLoop.size();
			patchF(i, 1) = i;
			patchF(i, 0) = (1 + i) % originalLoop.size();
		}

		// also upsample patch. patch and original mesh will have the same upsampling level now
		// making it trivial to fuse them together.
		igl::upsample(Eigen::MatrixXd(patchV), Eigen::MatrixXi(patchF), patchV, patchF, upsampleN);
	}

	// the final mesh, where the patch and original mesh has been fused together.
	std::vector<std::vector<double>> fusedV;
	std::vector<std::vector<int>> fusedF;

	int index = 0; // vertex index counter for the fused mesh.

				   // first, we add the upsampled patch to the fused mesh.
	{
		for (; index < patchV.rows(); ++index) {
			fusedV.push_back({ patchV(index, 0), patchV(index, 1), patchV(index, 2) });
		}

		int findex = 0;
		for (; findex < patchF.rows(); ++findex) {
			fusedF.push_back({ patchF(findex, 0), patchF(findex, 1), patchF(findex, 2) });
		}
	}

	// now, we fuse the patch and the original mesh together.
	{
		// remaps indices from the original mesh to the fused mesh.
		std::map<int, int> originalToFusedMap;

		for (int itri = 0; itri < originalF.rows(); ++itri) {

			int triIndices[3];
			for (int iv = 0; iv < 3; ++iv) {

				int triIndex = originalF(itri, iv);

				int ret;

				if (originalToFusedMap.count(triIndex) == 0) {
					int foundMatch = -1;

					// the vertices at the boundary are the same, for both the two meshes(patch and original mesh).
					// we ensure that these vertices are not duplicated.
					// this is also how we ensure that the two meshes are fused together.
					for (int jj = 0; jj < patchV.rows(); ++jj) {
						VectorXd u(3); u << fusedV[jj][0], fusedV[jj][1], fusedV[jj][2];
						VectorXd v(3); v << originalV(triIndex, 0), originalV(triIndex, 1), originalV(triIndex, 2);

						if ((u - v).norm() < 0.00001) {
							foundMatch = jj;
							break;
						}
					}

					if (foundMatch != -1) {
						originalToFusedMap[triIndex] = foundMatch;
						ret = foundMatch;
					}
					else {
						fusedV.push_back({ originalV(triIndex, 0), originalV(triIndex, 1), originalV(triIndex, 2) });
						originalToFusedMap[triIndex] = index;
						ret = index;
						index++;
					}
				}
				else {
					ret = originalToFusedMap[triIndex];
				}

				triIndices[iv] = ret;
			}

			fusedF.push_back({
				triIndices[0],
				triIndices[1],
				triIndices[2] });


		}

	}

	MatrixXd fairedV(fusedV.size(), 3);
	MatrixXi fairedF(fusedF.size(), 3);
	// now we shall do surface fairing on the mesh, to ensure
	// that the patch conforms to the surrounding curvature.
	{

		for (int vindex = 0; vindex < fusedV.size(); ++vindex) {
			auto r = fusedV[vindex];

			fairedV(vindex, 0) = r[0];
			fairedV(vindex, 1) = r[1];
			fairedV(vindex, 2) = r[2];
		}

		for (int findex = 0; findex < fusedF.size(); ++findex) {
			auto r = fusedF[findex];

			fairedF(findex, 0) = r[0];
			fairedF(findex, 1) = r[1];
			fairedF(findex, 2) = r[2];
		}

		VectorXi b(fairedV.rows() - patchV.rows());
		MatrixXd bc(fairedV.rows() - patchV.rows(), 3);
		// setup the boundary conditions. This is simply the vertex positions of the vertices not part of the patch.
		for (int i = (int)patchV.rows(); i < (int)fairedV.rows(); ++i) {
			int jj = i - (int)patchV.rows();

			b(jj) = i;

			bc(jj, 0) = fairedV(i, 0);
			bc(jj, 1) = fairedV(i, 1);
			bc(jj, 2) = fairedV(i, 2);
		}

		MatrixXd Z;
		int k = 2;
		// surface fairing simply means that we solve the equation
		// Delta^2 f = 0
		// with appropriate boundary conditions.
		// this function igl::harmonic from libigl takes care of that.

		// note that this is pretty inefficient thought.
		// the only boundary conditions necessary are the 2-ring of vertices around the patch.
		// the rest of the mesh vertices need not be specified.
		// we specify the rest of the mesh for simplicity of code, but it is not strictly necessary,
		// and pretty inefficient, since we have to solve a LARGE matrix equation as a result of this.
		igl::harmonic(fairedV, fairedF, b, bc, k, Z);
		fairedV = Z;
	}

	// finally, we do a decimation step.
	MatrixXd finalV(fusedV.size(), 3);
	MatrixXi finalF(fusedF.size(), 3);
	VectorXi temp0; VectorXi temp1;
	igl::decimate(fairedV, fairedF, outFacesN, finalV, finalF, temp0, temp1);
	
	igl::writeOFF(outFile, finalV, finalF);
}

