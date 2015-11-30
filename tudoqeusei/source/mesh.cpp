#pragma warning( disable : 4267)
#pragma comment(lib,"glew32.lib")
#include <GL/glew.h>
#include <GL/freeglut.h> 
#include <fstream>
#include "mesh.h"
#include "glsl_wrapper.h"
#include "math_headers.h"

/*------------------------------------------

Como essa classe funciona:  a classe Mesh é uma classe de malha genérica
							a classe Primitive é uma classe que deriva da Mesh

							Primitiva pode ser um cubo, esfera, ou um plano (só fiz essas primitivas xD)
							elas são geradas na hora, via equaçõeszinhas

							a classe Mesh serve mais é pra guardar modelos OBJ

Modelos 3D ficam armazenados em um vetor chamado m_positions
esse vetor é do tipo glm::vec3, entao cada posicao dele é um vetor de 3 dimensoes

de forma pratica, um modelo 3D é formado por N pontos com dimensao x,y e z
ligados de forma a criar triangulos.

como criar um ponto usando a glm

glm::vec3 p = glm::vec3(0,0,0) 
cria um ponto p na origem.

um modelo 3D é composto de seus vertices, e da triangulação
a triangulação é um vetor de inteiros que agem como indices no vetor das posições dos vertices.

Cada três inteiros formam um triangulo, que a função glDrawElements(GL_TRIANGLES, ...) usa

exemplo: quero desenhar um quadrado

um quadrado precisa de quatro pontos
glm::vec3 p1,p2,p3,p4;

p1 = glm::vec3(0,0,0);
p2 = glm::vec3(1,0,0);
p3 = glm::vec3(0,0,1);
p4 = glm::vec3(1,0,1);

//guardo no vetor de vertices
m_positions.push_back(p1);
m_positions.push_back(p2);
m_positions.push_back(p3);
m_positions.push_back(p4);

o minimo de triangulos para formar um quadrado são dois:

p3-----p4
|\     |
|  \   |
|    \ |
p1-----p2

bom, o primeiro triangulo é formado pelo ponto p1, p2, p3

m_triangles.push_back(0); <- aqui vai o indice de p1 no vetor m_positions
m_triangles.push_back(1);
m_triangles.push_back(2);

e o segundo triangulo, pelos pontos p3, p4 e p2

m_triangles.push_back(2);
m_triangles.push_back(3);
m_triangles.push_back(1);

ja que a conectividade nunca muda nos nossos modelos (afinal, nao estamos cortando eles)
podemos mudar tranquilamente a posição dos vertices e deformar ou mudar o modelo sem passarmos trabalho

exemplo:

m_positions[0] += glm::vec3(0,1,0);  -> desloquei p1 em 1 unidade na direção y (soma de vetores)


qualquer coisa ._. eu tento explicar mais...
----------------------------------------------*/

void Mesh::Reset()
{
    Cleanup();
    Init(); //to usando heranca, entao cada classe tem seu init
	computeNormal(); //computa as normais...
}
void Mesh::Cleanup()
{
	//limpa todo mundo
    m_edge_list.clear();
    m_positions.clear();
    m_normals.clear();
    m_colors.clear();
    m_triangle_list.clear();
	m_indices.clear();

	//Zera todo mundo
	m_current_positions.setZero();
	m_current_velocities.setZero();
	m_previous_velocities.setZero();

	m_mass_matrix.setZero();
	m_inv_mass_matrix.setZero();
}
void Mesh::Draw(const VBO& vbos, bool wire_frame, bool show_texture)
{
	//vamo la !

	//essa funcao computa a normal de cada face (triangulo) do modelo.
	//a normal é necessaria pra que o modelo seja visivel. 
	//entenda como se fosse o vetor que diz que o objeto "reflete" luz, e não absorve toda ela ^^
	computeNormal();

	//digo se o modelo vai ser só outline ou se é preenchido ^^
    glPolygonMode(GL_FRONT_AND_BACK, (wire_frame ? GL_LINE : GL_FILL));

	//numero de vertices
    unsigned int size = m_vertices_number;

	//numero de triangulos do modelo atual
    unsigned int element_num = m_triangle_list.size();

    //linko a posicao dos vertices pra placa de video
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);

	//se o objeto é  simulado fisicamente, colorizo ele com os vetores de força e atualizo as posições
	if (!m_mesh_type == MESH_TYPE_OBJ)
	{
		for (unsigned int i = 0; i < size; ++i)
		{
			m_positions[i] = glm::vec3(m_current_positions[3 * i + 0], m_current_positions[3 * i + 1], m_current_positions[3 * i + 2]);
		}
		for (unsigned int i = 0; i < size; ++i)
		{
			m_colors[i] = glm::vec3(m_current_velocities[3 * i + 0], m_current_velocities[3 * i + 1], m_current_velocities[3 * i + 2]);
		}
	}
    glBufferData(GL_ARRAY_BUFFER, 3 * size * sizeof(float), &m_positions[0], GL_DYNAMIC_DRAW);

    //linko a cor de cada vertice para a placa de video
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * size * sizeof(float), &m_colors[0], GL_STATIC_DRAW);
    
	//linko as normais de cada vertice para a placa de video
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
    glBufferData(GL_ARRAY_BUFFER, 3 * size * sizeof(float), &m_normals[0], GL_DYNAMIC_DRAW);
    
	//linko as coordenadas de textura de cada vertice para a placa de video
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_tbo);
    glBufferData(GL_ARRAY_BUFFER, 2 * size * sizeof(float), &m_texcoords[0], GL_STATIC_DRAW);

    //linko o indice dos vertices de cada triangulo para a placa de video
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, element_num * sizeof(unsigned int), &m_triangle_list[0], GL_STATIC_DRAW);

	//envio todo mundo 
    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_vbo);
	glEnableVertexAttribArray(RenderWrapper::VERTEX_COORD_ATTRIB);
    glVertexAttribPointer(RenderWrapper::VERTEX_COORD_ATTRIB, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_cbo);
	glEnableVertexAttribArray(RenderWrapper::COLOR_ATRIB);
    glVertexAttribPointer(RenderWrapper::COLOR_ATRIB, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_nbo);
	glEnableVertexAttribArray(RenderWrapper::NORMAL_ATTRIB);
    glVertexAttribPointer(RenderWrapper::NORMAL_ATTRIB, 3, GL_FLOAT, GL_FALSE, 0, 0);

    glBindBuffer(GL_ARRAY_BUFFER, vbos.m_tbo);
	glEnableVertexAttribArray(RenderWrapper::TEXTURE_COORD_ATTRIB);
	glVertexAttribPointer(RenderWrapper::TEXTURE_COORD_ATTRIB, 2, GL_FLOAT, GL_FALSE, 0, 0);

    glm::mat4 identity = glm::mat4(); // identity matrix
    glUniformMatrix4fv(vbos.m_uniform_transformation, 1, false, &identity[0][0]);

	//habilito ou desabilito a textura
    glUniform1i(vbos.m_uniform_enable_texture, show_texture);

	//desenho de fato os triangulos
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbos.m_ibo);
    glDrawElements(GL_TRIANGLES, element_num, GL_UNSIGNED_INT, 0);

	//desabilito tudo
    glDisableVertexAttribArray(RenderWrapper::VERTEX_COORD_ATTRIB);
    glDisableVertexAttribArray(RenderWrapper::COLOR_ATRIB);
    glDisableVertexAttribArray(RenderWrapper::NORMAL_ATTRIB);
    glDisableVertexAttribArray(RenderWrapper::TEXTURE_COORD_ATTRIB);

    glUniform1i(vbos.m_uniform_enable_texture, 0); // desabilita texture

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); 
}
void Mesh::Draw2()
{
	glBegin(GL_LINES);
	for (auto c : m_edge_list)
	{
		if (c.m_v1 == ((40 * 40) / 2) || c.m_v2 == ((40 * 40) / 2))
		{
			
				glColor3f(1.0, 1.0, 1.0);
				glVertex3f(m_positions[c.m_v1].x, m_positions[c.m_v1].y, m_positions[c.m_v1].z);
				glVertex3f(m_positions[c.m_v2].x, m_positions[c.m_v2].y, m_positions[c.m_v2].z);

		}

	}
	glEnd();
}
void Mesh::computeNormal()
{
	//explico isso via chamada se tu quiser...mas é algebra marota
    // reset all the normal.
    glm::vec3 zero(0.0);
    for(std::vector<glm::vec3>::iterator n = m_normals.begin(); n != m_normals.end(); ++n)
    {
        *n = zero;
    }
    // calculate normal for each individual triangle
    unsigned int triangle_num = m_triangle_list.size() / 3;
    unsigned int id0, id1, id2;
    glm::vec3 p0, p1, p2;
    glm::vec3 normal;
    for(unsigned int i = 0; i < triangle_num; ++i)
    {
        id0 = m_triangle_list[3 * i];
        id1 = m_triangle_list[3 * i + 1];
        id2 = m_triangle_list[3 * i + 2];

        p0 = m_positions[id0];
        p1 = m_positions[id1];
        p2 = m_positions[id2];

		glm::vec3 pa, pb;
		pa = p1-p0;
		pb = p2-p1;

        normal = glm::cross((pa),(pb));
		glm::vec3 normalN = normal;
		//normalN = glm::normalize(normalN);

        glm::vec3 glm_normal = glm::vec3(normalN[0], normalN[1], normalN[2]);

        m_normals[id0] += glm_normal;
        m_normals[id1] += glm_normal;
        m_normals[id2] += glm_normal;
    }
    // re-normalize all the normals.
    for(std::vector<glm::vec3>::iterator n = m_normals.begin(); n != m_normals.end(); ++n)
    {
        if (glm::length(*n) > 0.000001) // skip if norm is a zero vector
            *n = glm::normalize(*n);
    }

}
void Mesh::move_to(const glm::vec3& location)
{
	m_tdistance+=location;
	unsigned int index;
    for(index = 0; index < m_vertices_number; ++index)
    {
        m_positions[index] += location;
    }
}
void Mesh::scale_primitive(const float& factor)
{
	unsigned int index;
	for(index = 0; index < m_vertices_number; ++index)
	{
		m_positions[index] -= m_tdistance;
		m_positions[index] *= factor;
		m_positions[index] += m_tdistance;
	}
}
void Mesh::change_color(const glm::vec3& color)
{
    for (int i = 0; i < m_colors.size(); i++)
    {
        m_colors[i] = color;
    }
}


void Mesh::releaseMesh()
{
	mouseInteraction = false;
	grab_index = 0;
}

void Mesh::grabMesh(glm::vec3 pos)
{
	EigenVector3 position = GLM2Eigen(pos);
	ScalarType dist = (m_current_positions.block_vector(0) - position).norm();
	unsigned int current_index = 0;
	for (int i = 0; i < m_vertices_number; i++)
	{
		ScalarType dist2 = (m_current_positions.block_vector(i) - position).norm();
		if (dist2 < dist)
		{
			current_index = i;
			dist = dist2;

		}
	}
	grab_index = current_index;
}

void Mesh::calculatePointDisplacement(Camera c, glm::vec2 mouseOld, glm::vec2 mouseNew, float mouseZ)
{
	float z = (m_current_positions.block_vector(grab_index) - GLM2Eigen(c.m_position)).norm();
	float fatZ = (z)*0.0014f;

	mouse3DPosition = m_current_positions.block_vector(grab_index);
	mouse3DPosition += GLM2Eigen((mouseNew.y - mouseOld.y)*c.m_up)*-fatZ;
	mouse3DPosition += GLM2Eigen((mouseNew.x - mouseOld.x)*c.m_right)*fatZ;

	m_current_positions.block_vector(grab_index) = mouse3DPosition;
	m_current_velocities.block_vector(grab_index) = EigenVector3(0, 0, 0);
}

void Mesh::displacePoint()
{
	m_current_positions.block_vector(grab_index) = mouse3DPosition*0.1 + m_current_positions.block_vector(grab_index)*0.9;
	m_current_velocities.block_vector(grab_index) = EigenVector3(0, 0, 0);
}

//OBJ CLASS METHODS
bool ObjMesh::Init(){
	
    //m_loaded_mesh = new MeshLoader(m_tet_file_path, m_tet_scaling);
    if (read_from_file(filename) == false)
    {
        std::cout << "Error Loading OBJ." << std::endl;
        return false;
    }

	//initializing model
    generateParticleList();
    generateTriangleList();
	if (m_mesh_type == MESH_TYPE_OBJ_PHIS)
	{
		genericEdgeList();
		generateEdgeList();
	}

	//coelho é grande, então to escalando e movendo ele pra baixo em 0.3 unidades
	scale_primitive(0.1);
	move_to(glm::vec3(0,-0.3f,0));

	fprintf(stdout, "OBJ loaded...\n");
    return true;
}
bool ObjMesh::read_from_file(char* filename)
{
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();
    glm::vec3 mat_color(0.6); 

    std::ifstream infile(filename);
    if(!infile.good())
    {
        printf("Error in loading file %s\n", filename);
        return false;
    }
    char buffer[256];
    unsigned int ip0, ip1, ip2;
    unsigned int n0, n1, n2;
    glm::vec3 pos;
	float m_scaling = 1;
    while(!infile.getline(buffer,255).eof())
    {
        buffer[255] = '\0';
        if(buffer[0] == 'v' && (buffer[1] == ' ' || buffer[1] == 32))
        {
            if(sscanf_s(buffer, "v %f %f %f", &pos.x, &pos.y, &pos.z) == 3)
            {
                pos = m_scaling * pos;
                m_positions.push_back(pos);
            }
            else
            {
                printf("Vertex is not in desired format.\n");
                return false;
            }
        }
        else if (buffer[0] == 'v' && buffer[1] == 'n' && (buffer[2] == ' ' || buffer[2] == 32))
        {
            // load normals from obj file.
        }
        else if (buffer[0] == 'f' && (buffer[1] == ' ' || buffer[1] == 32))
        {
            if(sscanf_s(buffer, "f %u %u %u", &ip0, &ip1, &ip2) == 3)
            {
                m_indices.push_back(--ip0);
                m_indices.push_back(--ip1);
                m_indices.push_back(--ip2);
            }
            else if(sscanf_s(buffer, "f %u//%u %u//%u %u//%u", &ip0, &n0, &ip1, &n1, &ip2, &n2) == 6)
            {
                m_indices.push_back(--ip0);
                m_indices.push_back(--ip1);
                m_indices.push_back(--ip2);
            }
            else if(sscanf_s(buffer, "f %u/%u %u/%u %u/%u", &ip0, &n0, &ip1, &n1, &ip2, &n2) == 6)
            {
                m_indices.push_back(--ip0);
                m_indices.push_back(--ip1);
                m_indices.push_back(--ip2);
            }
            else
            {
                printf("Triangle indices is not in desired format.\n");
                return false;
            }
        }
    }
    // normals

    unsigned int id, size;
    bool vert_norm = (m_normals.size() != m_positions.size());
    if(vert_norm)
        m_normals.resize(m_positions.size(), glm::vec3(0.0f));

    size = m_indices.size();
    glm::uvec3 triangle;
    glm::vec3 p0, p1, p2;
    glm::vec3 norm;
    float phi0, phi1, phi2;
    float pi = glm::radians(180.0f);
    for(id = 0; id < size; id+=3)
    {
        triangle = glm::uvec3(m_indices[id], m_indices[id+1], m_indices[id+2]);
        p0 = m_positions[triangle.x];
        p1 = m_positions[triangle.y];
        p2 = m_positions[triangle.z];
        norm = glm::normalize(glm::cross(p1 - p0, p2 - p0));
        // calculate vertex normal
        if(vert_norm)
        {
            phi0 = std::acos(glm::dot(p1 - p0, p2 - p0) / (glm::length(p1 - p0) * glm::length(p2 - p0)));
            phi1 = std::acos(glm::dot(p0 - p1, p2 - p1) / (glm::length(p0 - p1) * glm::length(p2 - p1)));
            phi2 = pi - phi0 - phi1;

            m_normals[triangle.x] += phi0 * norm;
            m_normals[triangle.y] += phi1 * norm;
            m_normals[triangle.z] += phi2 * norm;
        }
    }
    // re-normalize all normals
    for(std::vector<glm::vec3>::iterator iter = m_normals.begin(); iter != m_normals.end(); ++iter)
    {
        *iter = glm::normalize(*iter);
        m_colors.push_back(mat_color);
    }
	return true;
}
bool Mesh::write_to_file(char* filename)
{
	std::ofstream outfile;
	outfile.open(filename, std::ifstream::out);
	if (outfile.is_open())
	{
		outfile<<"#Wgrampon OBJ exporter"<<std::endl;
		glm::vec3 pos, tri_vert;
		for(int i =0; i<m_vertices_number; i++)
		{
			pos = m_positions[i];
			outfile<<"v "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
		}
		outfile<<"# "<<m_vertices_number<<" vertices"<<std::endl;

		for(int i=0; i<m_triangle_list.size()/3; i++)
		{
			tri_vert = glm::vec3(m_triangle_list[3*i]+1,m_triangle_list[3*i+1]+1,m_triangle_list[3*i+2]+1);
			outfile<<"f "<<tri_vert[0]<<" "<<tri_vert[1]<<" "<<tri_vert[2]<<std::endl;
		}
		outfile<<"# "<<m_triangle_list.size()/3<<" triangles"<<std::endl;
		outfile.close();
		return true;
	}
	else
	{
		std::cout<<"Error writting obj"<<std::endl;
		return false;
	}
}
void ObjMesh::generateParticleList()
{
    m_vertices_number = m_positions.size();
    m_system_dimension = m_vertices_number * 3;
    m_positions.resize(m_vertices_number);
    m_normals.resize(m_vertices_number);
    m_colors.resize(m_vertices_number);

    // color
    change_color(glm::vec3(0.47,0.08,0.07));

	m_total_mass = 1;
	ScalarType unit_mass = m_total_mass / m_vertices_number;

	m_current_positions.resize(m_system_dimension);
	m_current_velocities.resize(m_system_dimension);
	m_previous_velocities.resize(m_system_dimension);

	m_mass_matrix.resize(m_system_dimension, m_system_dimension);
	m_inv_mass_matrix.resize(m_system_dimension, m_system_dimension);
	m_identity_matrix.resize(m_system_dimension, m_system_dimension);

	m_current_positions.setZero();

	unsigned int index;
	for (index = 0; index < m_vertices_number; ++index)
	{
		m_current_positions.block_vector(index) = GLM2Eigen(m_positions[index]);
	}
	// Assign initial velocity to zero
	m_current_velocities.setZero();
	m_previous_velocities.setZero();

	// Assign mass matrix and an equally sized identity matrix
	std::vector<SparseMatrixTriplet> i_triplets;
	std::vector<SparseMatrixTriplet> m_triplets;
	std::vector<SparseMatrixTriplet> m_inv_triplets;
	i_triplets.clear();
	m_triplets.clear();
	ScalarType inv_unit_mass = 1.0 / unit_mass;
	for (index = 0; index < m_system_dimension; index++)
	{
		i_triplets.push_back(SparseMatrixTriplet(index, index, 1));
		m_triplets.push_back(SparseMatrixTriplet(index, index, unit_mass));
		m_inv_triplets.push_back(SparseMatrixTriplet(index, index, inv_unit_mass));
	}
	m_mass_matrix.setFromTriplets(m_triplets.begin(), m_triplets.end());
	m_inv_mass_matrix.setFromTriplets(m_inv_triplets.begin(), m_inv_triplets.end());
	m_identity_matrix.setFromTriplets(i_triplets.begin(), i_triplets.end());

}
void ObjMesh::generateTriangleList()
{
	m_triangle_list.resize(m_indices.size());
	for (int i = 0; i != m_indices.size(); ++i)
    {
        m_triangle_list[i] = m_indices[i];
    }
}

//void Mesh::generateEdgeList()
//{
//	// generate all the edges from the vertices and triangle list.
//	// courtesy of Eric Lengyel, "Building an Edge List for an Arbitrary Mesh". Terathon Software 3D Graphics Library.
//	// http://www.terathon.com/code/edges.html
//	unsigned int vert_num = m_vertices_number;
//	unsigned int tri_num = m_triangle_list.size() / 3;
//
//	unsigned int *first_edge = new unsigned int[vert_num + 3 * tri_num];
//	unsigned int *next_edge = first_edge + vert_num;
//
//	for (unsigned int i = 0; i < vert_num; ++i)
//		first_edge[i] = 0xFFFFFFFF;
//	// First pass over all triangles. Finds out all the edges satisfying the condition that
//	// the first vertex index is less than the second vertex index when the direction from 
//	// the first to the second represents a counterclockwise winding around the triangle to
//	// which the edge belongs. For each edge found, the edge index is stored in a linked 
//	// list of edges belonging to the lower-numbered vertex index i. This allows us to 
//	// quickly find an edge in the second pass whose higher-numbered vertex is i.
//
//	unsigned int edge_count = 0;
//	const unsigned int* triangle = &m_triangle_list[0];
//	unsigned int i1, i2;
//	for (unsigned int t = 0; t < tri_num; ++t)
//	{
//		i1 = triangle[2];
//		for (unsigned int n = 0; n < 3; ++n)
//		{
//			i2 = triangle[n];
//			if (i1 < i2)
//			{
//				Edge new_edge;
//				new_edge.m_v1 = i1;
//				new_edge.m_v2 = i2;
//				new_edge.m_tri1 = t;
//				new_edge.m_tri2 = t;
//				m_edge_list.push_back(new_edge);
//
//				unsigned int edge_idx = first_edge[i1];
//				if (edge_idx == 0xFFFFFFFF)
//				{
//					first_edge[i1] = edge_count;
//				}
//				else
//				{
//					while (true)
//					{
//						unsigned int idx = next_edge[edge_idx];
//						if (idx == 0xFFFFFFFF)
//						{
//							next_edge[edge_idx] = edge_count;
//							break;
//						}
//						edge_idx = idx;
//					}
//				}
//
//				next_edge[edge_count] = 0xFFFFFFFF;
//				edge_count++;
//			}
//			i1 = i2;
//		}
//		triangle += 3;
//	}
//
//	// Second pass over all triangles. Finds out all the edges satisfying the condition that
//	// the first vertex index is greater than the second vertex index when the direction from 
//	// the first to the second represents a counterclockwise winding around the triangle to
//	// which the edge belongs. For each of these edges, the same edge should have already been
//	// found in the first pass for a different triangle. So we search the list of edges for the
//	// higher-numbered index for the matching edge and fill in the second triangle index. The 
//	// maximum number of the comparisons in this search for any vertex is the number of edges
//	// having that vertex as an endpoint.
//	triangle = &m_triangle_list[0];
//	for (unsigned int t = 0; t < tri_num; ++t)
//	{
//		i1 = triangle[2];
//		for (unsigned int n = 0; n < 3; ++n)
//		{
//			i2 = triangle[n];
//			if (i1 > i2)
//			{
//				bool is_new_edge = true;
//				for (unsigned int edge_idx = first_edge[i2]; edge_idx != 0xFFFFFFFF; edge_idx = next_edge[edge_idx])
//				{
//					Edge *edge = &m_edge_list[edge_idx];
//					if ((edge->m_v2 == i1) && (edge->m_tri1 == edge->m_tri2))
//					{
//						edge->m_tri2 = t;
//						is_new_edge = false;
//						break;
//					}
//				}
//				// for case where a edge belongs to only one triangle. i.e. mesh is not watertight.
//				if (is_new_edge)
//				{
//					Edge new_edge;
//					new_edge.m_v1 = i1;
//					new_edge.m_v2 = i2;
//					new_edge.m_tri1 = t;
//					new_edge.m_tri2 = t;
//					m_edge_list.push_back(new_edge);
//
//					unsigned int edge_idx = first_edge[i1];
//					if (edge_idx == 0xFFFFFFFF)
//					{
//						first_edge[i1] = edge_count;
//					}
//					else
//					{
//						while (true)
//						{
//							unsigned int idx = next_edge[edge_idx];
//							if (idx == 0xFFFFFFFF)
//							{
//								next_edge[edge_idx] = edge_count;
//								break;
//							}
//							edge_idx = idx;
//						}
//					}
//
//					next_edge[edge_count] = 0xFFFFFFFF;
//					edge_count++;
//				}
//			}
//			i1 = i2;
//		}
//		triangle += 3;
//	}
//
//	delete[] first_edge;
//	printf("Edge number: %u.\n", m_edge_list.size());
//}


//phisics related ?
void Mesh::genericEdgeList()
{
	for (int i = 0; i < m_triangle_list.size() / 3; i++)
	{
		Edge e1, e2, e3;
		e1.m_v1 = m_triangle_list[3 * i];
		e1.m_v2 = m_triangle_list[3 * i + 1];
		
		e2.m_v1 = m_triangle_list[3 * i];
		e2.m_v2 = m_triangle_list[3 * i + 2];

		e3.m_v1 = m_triangle_list[3 * i + 1];
		e3.m_v2 = m_triangle_list[3 * i + 2];

		m_edge_list_temp.push_back(e1);
		m_edge_list_temp.push_back(e2);
		m_edge_list_temp.push_back(e3);
	}
}
void Mesh::generateEdgeList()
{
	unsigned int i1, i2;
	unsigned int vert_num = m_positions.size();

	SparseMatrix EdgeMatrix(vert_num, vert_num);
	EdgeMatrix.setZero();

	for (int i = 0; i<m_edge_list_temp.size(); ++i)
	{
		i1 = m_edge_list_temp[i].m_v1;
		i2 = m_edge_list_temp[i].m_v2;

		if (EdgeMatrix.coeff(i1, i2) < EPSILON)
		{
			EdgeMatrix.coeffRef(i1, i2) = 1;
			EdgeMatrix.coeffRef(i2, i1) = 1;

			Edge new_edge;
			new_edge.m_v1 = i1;
			new_edge.m_v2 = i2;
			new_edge.l = glm::length(m_positions[i2] - m_positions[i1]);
			m_edge_list.push_back(new_edge);
		}
	}
	m_edge_list_temp.clear();
	std::cout << "Edge list complete" << std::endl;
}