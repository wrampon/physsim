#include <fstream>
#include "primitive.h"
#include <time.h>

#define COLLISION_EPSILON 1e-3

//----------Base Class-----------//
void Primitive::Draw2(){
	glBegin(GL_LINES);
	for (auto c : m_edge_list)
	{
		if (c.m_v1 == ((40 * 40) / 2)+20 || c.m_v2 == ((40 * 40) / 2)+20)
		{

			glColor3f(1.0, 1.0, 1.0);
			glVertex3f(m_positions[c.m_v1].x, m_positions[c.m_v1].y, m_positions[c.m_v1].z);
			glVertex3f(m_positions[c.m_v2].x, m_positions[c.m_v2].y, m_positions[c.m_v2].z);

		}

	}
	glEnd();
}
void Primitive::generateTriangleList()
{
	m_triangle_list.resize(m_indices.size());
	for (int i = 0; i != m_indices.size(); ++i)
    {
        m_triangle_list[i] = m_indices[i];
    }

}

//----------Plane Class----------//
bool Plane::Init()
{
	srand (time(NULL));
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    glm::vec3 mat_color(0.6);
	double norm = 1;
	double increment = 2*norm/(side);
	glm::vec3 point;
	unsigned int i = 0;
	double j,k,l, it1, it2;
	j = k = l = -norm;
	int sideAm = 0;

	int upLimit = side;
	double lowerlimitY = -norm+increment;
	double lowerlimitX = lowerlimitY;

	for( it1 = 0 , k = lowerlimitY; it1<upLimit; k+=increment, it1++)
	{
		sideAm++;
		for( it2 = 0, l = lowerlimitX; it2<upLimit; l+=increment, it2++)
		{
			point = glm::vec3(l,-1,k);
			m_positions.push_back(point);
			m_colors.push_back(mat_color);
		}
	}

	side = sideAm;

	//integration variables
	m_vertices_number =  m_positions.size();
	m_system_dimension =  m_vertices_number * 3; 
	std::cout << m_system_dimension << std::endl;
	 // Assign initial position to all the vertices.
	m_current_positions.resize(m_system_dimension);
	m_current_velocities.resize(m_system_dimension);
	m_previous_velocities.resize(m_system_dimension);

	m_mass_matrix.resize(m_system_dimension, m_system_dimension);
	m_inv_mass_matrix.resize(m_system_dimension, m_system_dimension);
	m_identity_matrix.resize(m_system_dimension, m_system_dimension);

	m_current_positions.setZero();
	unsigned int b, c, index;
	for (b = 0; b < m_positions.size(); ++b)
	{
		m_current_positions.block_vector(b) = EigenVector3(m_positions[b].x, m_positions[b].y, m_positions[b].z);
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

	m_total_mass = 1;
	ScalarType unit_mass = m_total_mass / m_vertices_number;
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

	generateTriangleList();
	generateEdgeList();

	scale_primitive(0.5);
	m_normals.resize(m_positions.size());
	change_color(glm::vec3(0.9,0.9,0.9));

	float coord = (float)1.0/(float)(side-1);
	//texture assignment
	for(int i = 0; i<side; i++)
		for(int j = 0; j<side; j++)
			m_texcoords.push_back(glm::vec2(j*coord, i*coord));
	
	return true;
}

void Plane::generateTriangleList()
{
	m_indices.clear();
	for(unsigned int  k = 0; k<(side-1); k++)
	{
		for(unsigned int l = 0; l<(side-1); l++)
		{
			//change normal to outside according to face of plane
			m_indices.push_back(k*side + l);
			m_indices.push_back((k+1)*side + l);
			m_indices.push_back(k*side + l+1);

			m_indices.push_back((k+1)*side + l+1);
			m_indices.push_back(k*side + l+1);
			m_indices.push_back((k+1)*side + l);

			Edge e1, e2, e3, e4, e5, e6, e7, e8;
			//up
			e1.m_v1 = k*side + l;
			e1.m_v2 = (k + 1)*side + l;
			e1.t = STRETCH;

			//right
			e2.m_v1 = k*side + l;
			e2.m_v2 = k*side + l + 1;
			e2.t = STRETCH;

			//diagonal right
			e3.m_v1 = k*side + l;
			e3.m_v2 = (k + 1)*side + l +1;
			e3.t = STRETCH;

			//diagonal left
			if (l > 0)
			{
				e4.m_v1 = k*side + l;
				e4.m_v2 = (k + 1)*side + l - 1;
				e4.t = STRETCH;
				m_edge_list_temp.push_back(e4);

				//diagonal left + 2
				if (l > 1 && k<(side - 2))
				{
					e5.m_v1 = k*side + l;
					e5.m_v2 = (k + 2)*side + l - 2;
					e5.t = STRETCH;
					m_edge_list_temp.push_back(e5);

				}
			}
			//right + 2
			if (l < (side - 2))
			{
				e6.m_v1 = k*side + l;
				e6.m_v2 = k*side + l + 2;
				e6.t = STRETCH;
				m_edge_list_temp.push_back(e6);
			}

			if (k < (side - 2))
			{
				e7.m_v1 = k*side + l;
				e7.m_v2 = (k+2)*side + l;
				e7.t = STRETCH;
				m_edge_list_temp.push_back(e7);
			}

			if (k < (side - 2) && l<(side-2))
			{
				e8.m_v1 = k*side + l;
				e8.m_v2 = (k + 2)*side + l+2;
				e8.t = STRETCH;
				m_edge_list_temp.push_back(e8);
			}


			m_edge_list_temp.push_back(e1);
			m_edge_list_temp.push_back(e2);
			m_edge_list_temp.push_back(e3);


		}
	}
	m_triangle_list.resize(m_indices.size());
	for (int i = 0; i < m_indices.size(); ++i)
    {
        m_triangle_list[i] = m_indices[i];
    }	
}
void Plane::createMountains()
{
	//metodo pra deformar o modelo...
	//no caso, serve mais é pro plano
	int index = rand()%m_positions.size();
	float radius = (float)rand()/(float)(RAND_MAX);
	radius*=0.25f;

	glm::vec3 center = m_positions[index];
	glm::vec3 deformation;
	glm::vec3 partial;
	for(int i = 0; i<m_positions.size();i++)
	{
		partial = m_positions[i] - center;
		deformation.y = (radius*radius) - ((partial.x*partial.x) + (partial.z*partial.z));
		if(deformation.y>0) 
		{
			m_positions[i] += glm::vec3(0,deformation.y,0)*0.7f;
		}
	}

}
//----------Sphere Class----------//
bool Sphere::Init()
{
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    glm::vec3 mat_color(0.6);
    unsigned int slice = 48, stack = 48;

    glm::vec3 tnormal(0.0, 1.0, 0.0), tpos;
    tpos = m_radius * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

    float theta_z, theta_y, sin_z;
    float delta_y = 360.0 / slice, delta_z = 180.0 / stack;
    //loop over the sphere
    for(theta_z = delta_z; theta_z < 179.99; theta_z += delta_z)
    {
        for(theta_y = 0.0; theta_y < 359.99; theta_y += delta_y)
        {
            sin_z = sin(glm::radians(theta_z));
            
            tnormal.x = sin_z * cos(glm::radians(theta_y));
            tnormal.y = cos(glm::radians(theta_z));
            tnormal.z = -sin_z * sin(glm::radians(theta_y));

            tpos = m_radius * tnormal;

            m_positions.push_back(tpos);
            m_normals.push_back(tnormal);
            m_colors.push_back(mat_color);
        }
    }
    tnormal = glm::vec3(0.0, -1.0, 0.0);
    tpos = m_radius * tnormal;

    m_positions.push_back(tpos);
    m_normals.push_back(tnormal);
    m_colors.push_back(mat_color);

    //indices
    unsigned int j = 0, k = 0;
    for(j = 0; j < slice - 1; ++j)
    {
        m_indices.push_back(0);
        m_indices.push_back(j + 1);
        m_indices.push_back(j + 2);
    }
    m_indices.push_back(0);
    m_indices.push_back(slice);
    m_indices.push_back(1);

    for(j = 0; j < stack - 2; ++j)
    {
        for(k = 1 + slice * j; k < slice * (j + 1); ++k)
        {
            m_indices.push_back(k);
            m_indices.push_back(k + slice);
            m_indices.push_back(k + slice + 1);

            m_indices.push_back(k);
            m_indices.push_back(k + slice + 1);
            m_indices.push_back(k + 1);
        }
        m_indices.push_back(k);
        m_indices.push_back(k + slice);
        m_indices.push_back(k + 1);

        m_indices.push_back(k);
        m_indices.push_back(k + 1);
        m_indices.push_back(k + 1 - slice);
    }

    unsigned int bottom_id = (stack - 1) * slice + 1;
    unsigned int offset = bottom_id - slice;
    for(j = 0; j < slice - 1; ++j)
    {
        m_indices.push_back(j + offset);
        m_indices.push_back(bottom_id);
        m_indices.push_back(j + offset + 1);
    }
    m_indices.push_back(bottom_id - 1);
    m_indices.push_back(bottom_id);
    m_indices.push_back(offset);

    if(m_indices.size() != 6 * (stack - 1) * slice)
        printf("indices number not correct!\n");
	
	m_vertices_number =  m_positions.size();
	m_system_dimension =  m_vertices_number * 3; 
	
	generateEdgeList();
	generateTriangleList();
	change_color(glm::vec3(0.17,0.02,0.01));
	return true;
}
void Sphere::generateEdgeList()
{
	std::cout<<"Edge list complete"<<std::endl;
}

//----------Cube Class----------//
bool Cube::Init()
{
    m_positions.clear();
    m_colors.clear();
    m_normals.clear();
    m_indices.clear();

    glm::vec3 mat_color(0.6);
	double norm = 1;
	double increment = 2*norm/(side);
	glm::vec3 point;
	unsigned int i = 0;
	double j,k,l, it1, it2;
	j = k = l = -norm;
	int sideAm = 0;

	int upLimit = side - 2;
	double lowerlimitY = -norm+increment;
	double lowerlimitX = lowerlimitY;
	for(int z = 0; z<6; z++)
	{

		for( it1 = 0 , k = lowerlimitY; it1<=upLimit; k+=increment, it1++)
		{
			if(z==0)sideAm++;
			for( it2 = 0, l = lowerlimitX; it2<=upLimit; l+=increment, it2++){
				switch (z)
				{
					case 0: point = glm::vec3(l,k,-norm); break;
					case 1: point = glm::vec3(l,k,+norm); break;
					case 2: point = glm::vec3(-norm,k,l); break;
					case 3: point = glm::vec3(+norm,k,l); break;
					case 4: point = glm::vec3(l,-norm,k); break;
					case 5: point = glm::vec3(l,+norm,k); break;
				}
				m_positions.push_back(point);
				m_colors.push_back(mat_color);
			}
		}
	}

	//integration variables
	m_vertices_number =  m_positions.size();
	m_system_dimension =  m_vertices_number * 3; 
	
	 // Assign initial position to all the vertices.
	generateEdgeList();
	generateTriangleList();
	m_normals.resize(m_positions.size());
	scale_primitive(0.05f);
	change_color(glm::vec3(0.47,0.08,0.07));
	return true;
}
void Cube::generateEdgeList()
{
	std::cout<<"Edge list complete"<<std::endl;
}
void Cube::generateTriangleList()
{
	for(unsigned int z = 0; z<=5; z++)
		for(unsigned int  k = 0; k<(side-1); k++)
			for(unsigned int l = 0; l<(side-1); l++)
			{
				//change normal to outside according to face of cube
				if(z==0 || z==3 || z==5)
				{
					m_indices.push_back(z*side*side + k*side + l);
					m_indices.push_back(z*side*side + (k+1)*side + l);
					m_indices.push_back(z*side*side + k*side + l+1);

					m_indices.push_back(z*side*side + (k+1)*side + l+1);
					m_indices.push_back(z*side*side + k*side + l+1);
					m_indices.push_back(z*side*side + (k+1)*side + l);
				}
				else
				{
					m_indices.push_back(z*side*side + k*side + l);													
					m_indices.push_back(z*side*side + k*side + l+1);		
					m_indices.push_back(z*side*side + (k+1)*side + l);	

					m_indices.push_back(z*side*side + (k+1)*side + l+1);	
					m_indices.push_back(z*side*side + (k+1)*side + l);	
					m_indices.push_back(z*side*side + k*side + l+1);	
				}

			}

	for(unsigned int l = 0; l<side-1; l++)
	{
		//bottom
		m_indices.push_back(l);
		m_indices.push_back(4*side*side +l+1);
		m_indices.push_back(4*side*side +l);		

		m_indices.push_back(4*side*side +l+1);
		m_indices.push_back(l);
		m_indices.push_back(l+1);

		//up
		m_indices.push_back(side*(side-1) + l);
		m_indices.push_back(5*side*side +l);
		m_indices.push_back(5*side*side +l+1);

		m_indices.push_back(5*side*side +l+1);
		m_indices.push_back(side*(side-1) + l+1);
		m_indices.push_back(side*(side-1) + l);

		//left
		m_indices.push_back(side*l);													
		m_indices.push_back(2*side*side +side*l);		
		m_indices.push_back(2*side*side +side*(l+1));

		m_indices.push_back(2*side*side +side*(l+1));		
		m_indices.push_back(side*(l+1));			
		m_indices.push_back(side*l);			
				
		//right
		m_indices.push_back(side*l + (side-1));
		m_indices.push_back(3*side*side +side*(l+1));
		m_indices.push_back(3*side*side +side*l);

		m_indices.push_back(3*side*side +side*(l+1));		
		m_indices.push_back(side*l + (side-1));
		m_indices.push_back(side*(l+1) + (side-1));

		//backface
		m_indices.push_back(side*side + l);													
		m_indices.push_back(4*side*side +side*(side-1) + l);		
		m_indices.push_back(4*side*side +side*(side-1) + (l+1));		

		m_indices.push_back(4*side*side +side*(side-1) + (l+1));		
		m_indices.push_back(side*side + l+1);								
		m_indices.push_back(side*side + l);

		m_indices.push_back(side*side +side*(side-1) + l);
		m_indices.push_back(5*side*side +side*(side-1) + l+1);
		m_indices.push_back(5*side*side +side*(side-1) + l);

		m_indices.push_back(5*side*side +side*(side-1) + l+1);
		m_indices.push_back(side*side +side*(side-1) + l);
		m_indices.push_back(side*side +side*(side-1) + l+1);
	
		//back left right
		m_indices.push_back(side*side + side*l);
		m_indices.push_back(2*side*side +side*(l+1)+side-1);
		m_indices.push_back(2*side*side +side*l+side-1);

		m_indices.push_back(2*side*side +side*(l+1)+side-1);
		m_indices.push_back(side*side + side*l);
		m_indices.push_back(side*side + side*(l+1));

		m_indices.push_back(side*side + side*l + (side-1));
		m_indices.push_back(3*side*side +side*l+side-1);
		m_indices.push_back(3*side*side +side*(l+1)+side-1);

		m_indices.push_back(3*side*side +side*(l+1)+side-1);
		m_indices.push_back(side*side + side*(l+1) + (side-1));
		m_indices.push_back(side*side + side*l + (side-1));

		////left face
		m_indices.push_back(2*side*side +side*(side-1) + l);
		m_indices.push_back(5*side*side +side*(l+1));
		m_indices.push_back(5*side*side +side*l);

		m_indices.push_back(5*side*side +side*(l+1));
		m_indices.push_back(2*side*side +side*(side-1) + l);
		m_indices.push_back(2*side*side +side*(side-1) + l+1);

		m_indices.push_back(2*side*side + l);
		m_indices.push_back(4*side*side + side*l);
		m_indices.push_back(4*side*side + side*(l+1));

		m_indices.push_back(4*side*side + side*(l+1));
		m_indices.push_back(2*side*side + l+1);
		m_indices.push_back(2*side*side + l);

		//right face
		m_indices.push_back(3*side*side +side*(side-1) + l);
		m_indices.push_back(5*side*side +side*l +(side-1));
		m_indices.push_back(5*side*side +side*(l+1) +(side-1));

		m_indices.push_back(5*side*side +side*(l+1) +(side-1));
		m_indices.push_back(3*side*side +side*(side-1) + l+1);
		m_indices.push_back(3*side*side +side*(side-1) + l);

		m_indices.push_back(3*side*side + l);
		m_indices.push_back(4*side*side + side*(l+1) +(side-1));
		m_indices.push_back(4*side*side + side*l +(side-1));

		m_indices.push_back(4*side*side + side*(l+1) +(side-1));
		m_indices.push_back(3*side*side + l);
		m_indices.push_back(3*side*side + l+1);
				
	}

	//extreme corners
	//back extreme
	m_indices.push_back(0);
	m_indices.push_back(4*side*side);
	m_indices.push_back(2*side*side);

	m_indices.push_back(side-1);
	m_indices.push_back(3*side*side);
	m_indices.push_back(4*side*side + side-1);

	m_indices.push_back((side-1)*side);
	m_indices.push_back(2*side*side + side*(side-1));
	m_indices.push_back(5*side*side);

	m_indices.push_back((side-1)*side + side-1);
	m_indices.push_back(5*side*side + side-1);
	m_indices.push_back(3*side*side + side*(side-1));

	//front
	m_indices.push_back(side*side);
	m_indices.push_back(2*side*side + side-1);
	m_indices.push_back(4*side*side + (side-1)*side);

	m_indices.push_back(side*side + side-1);
	m_indices.push_back(4*side*side + side*(side-1) + side-1);
	m_indices.push_back(3*side*side + side-1);

	m_indices.push_back(side*side + (side-1)*side);
	m_indices.push_back(5*side*side + side*(side-1));
	m_indices.push_back(2*side*side + side*(side-1) + side-1);

	m_indices.push_back(side*side + (side-1)*side + side-1);
	m_indices.push_back(3*side*side + side*(side-1) + side-1);
	m_indices.push_back(5*side*side + side*(side-1) + side-1);

	m_triangle_list.resize(m_indices.size());
	for (int i = 0; i != m_indices.size(); ++i)
    {
        m_triangle_list[i] = m_indices[i];
    }

}