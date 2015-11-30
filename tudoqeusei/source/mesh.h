#ifndef _MESH_H_
#define _MESH_H_

#include <vector>

#include "opengl_headers.h"
#include "math_headers.h"
#include "camera.h"

//forward declarations
class Camera;
class Simulation;

typedef enum
{
	MESH_TYPE_OBJ,
	MESH_TYPE_OBJ_PHIS,
	MESH_TYPE_PRIM,
	MESH_TYPE_CLOTH
} MeshType;

typedef enum
{
	STRETCH,
	BENDING
} EdgeType;

struct Edge
{
    unsigned int m_v1, m_v2; // indices of endpoint vertices
    unsigned int m_tri1, m_tri2; // indices of adjacent faces
	float l;
	EdgeType t;
};

struct Triangle
{
	int index_v1, index_v2, index_v3;
};

struct TriangleList
{
	std::vector<Triangle> m_triangle_list;
};

class Mesh
{
    friend class Camera;
	friend class Simulation;

public:
	VectorX m_current_positions; // 1x3m

    Mesh() : m_mesh_type() {}
    Mesh(MeshType mesh_type) : m_mesh_type(mesh_type) {}
    virtual ~Mesh() {Cleanup();}

    void Reset();
    virtual bool Init() {std::cout << "Warning: reach base class virtual init function." << std::endl; return false;}
    virtual void Cleanup();
	bool write_to_file(char* filename);

    virtual void Draw(const VBO& vbos, bool wire_frame = false, bool g_show_texture = false);
	virtual void Draw2();

	virtual void move_to(const glm::vec3& location);
	virtual void scale_primitive(const float& factor);
	virtual void change_color(const glm::vec3& color);
	virtual void calculatePointDisplacement(Camera c, glm::vec2 mouseOld, glm::vec2 mouseNew, float mouseZ);
	virtual void displacePoint();

	virtual void grabMesh(glm::vec3 position);
	virtual void releaseMesh();
	MeshType m_mesh_type;

protected:
   

    unsigned int m_vertices_number; // m
    unsigned int m_system_dimension; // 3m
	unsigned int grab_index;

   	glm::vec3 m_tdistance;
	ScalarType m_total_mass;

    std::vector<Edge> m_edge_list;
	std::vector<Edge> m_edge_list_temp;

	std::vector<glm::vec3> m_positions;
    std::vector<glm::vec3> m_normals;
    std::vector<glm::vec3> m_colors;
	std::vector<glm::vec2> m_texcoords;

    std::vector<unsigned int> m_triangle_list;
	std::vector<TriangleList> m_vertex_triangles;

	//for obj
	std::vector<unsigned int> m_indices;

	//for cloth
	unsigned int m_dim[2]; // width and length

	// vertices positions/previous positions/mass
	
	VectorX m_current_velocities;		// 1x3m
	VectorX m_previous_velocities;		// 1x3m
	SparseMatrix m_mass_matrix;			// 3mx3m
	SparseMatrix m_inv_mass_matrix;		// 3mx3m
	SparseMatrix m_identity_matrix;		// 3mx3m

	//mouse
	EigenVector3 mouse3DPosition;
	bool mouseInteraction;

protected:
    // initialize every particle pos / vel / mass / color.
    virtual void generateParticleList() {std::cout << "Warning: reach base class virtual function." << std::endl;}
    // generate triangle list from vetices
    virtual void generateTriangleList() {std::cout << "Warning: reach base class virtual function." << std::endl;}
    // generate edge list from the geometry representation.
	virtual void generateEdgeList();
	
	//desnecessario porenquanto	
	virtual void genericEdgeList();

    void computeNormal();
};

class ObjMesh : public Mesh
{

public:
	char* filename;
	ObjMesh(char* filename) : Mesh(MESH_TYPE_OBJ), filename(filename) {}
	virtual ~ObjMesh() {}
	virtual bool Init();
	
protected:

protected:
	bool read_from_file(char* filename);
	// initialize every particle pos / vel / mass / color.
	virtual void generateParticleList();
	// generate triangle list from vetices
	virtual void generateTriangleList();

};

class ClothMesh : public Mesh
{

public:
	ClothMesh() : Mesh(MESH_TYPE_CLOTH) {}
	ClothMesh(unsigned int dim0, unsigned int dim1) : Mesh(MESH_TYPE_CLOTH) { m_dim[0] = dim0; m_dim[1] = dim1; }
	virtual ~ClothMesh() {}

	virtual bool Init();

protected:

	// initialize every particle pos / vel / mass / color.
	virtual void generateParticleList();
	// generate triangle list from vetices
	virtual void generateTriangleList();
	// generate edge list from the geometry representation.
	virtual void generateEdgeList();
};

#endif
