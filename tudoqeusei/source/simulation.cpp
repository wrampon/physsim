#pragma warning( disable : 4996)

#include "simulation.h"

Simulation::Simulation()

{
}

Simulation::~Simulation()
{
    clearConstraints();
}

void Simulation::Reset()
{    
    m_inertia_y.resize(m_mesh->m_system_dimension);
    m_external_force.resize(m_mesh->m_system_dimension);

    setupConstraints();
    m_selected_attachment_constraint = NULL;
	m_h = 0.03;
	m_gravity_constant = 10;
	m_stiffness_stretch = 50;
	m_damping_coefficient = 0.2;
	m_resistance_coefficient = 0.1;
	g_gravity = true;
}

void Simulation::Update()
{

    // update external force
	if (g_gravity) calculateExternalForce();
	switch (m_integration_method)
	{
		case INTEGRATION_EXPLICIT_EULER_COMMON:
			integrateClassicExplicitEuler();
			break;
		case INTEGRATION_EXPLICIT_EULER_DISCRETE:
			integrateExplicitEuler();
			break;
		case INTEGRATION_EXPLICIT_EULER_SYMPLETIC:
			integrateExplicitEuler();
			break;
	}
    // update velocity and damp
}


void Simulation::DrawConstraints(const VBO& vbos)
{
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->Draw(vbos);
    }
}


void Simulation::clearConstraints()
{
    for (unsigned int i = 0; i < m_constraints.size(); ++i)
    {
        delete m_constraints[i];
    }
    m_constraints.clear();
}

void Simulation::setupConstraints()
{
    clearConstraints();

    switch(m_mesh->m_mesh_type)
    {
		case MESH_TYPE_OBJ_PHIS:
		{
			// generate stretch constraints. assign a stretch constraint for each edge.
			EigenVector3 p1, p2;
			for(std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
			{
				p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
				p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
				SpringConstraint *c = new SpringConstraint(&m_stiffness_stretch, e->m_v1, e->m_v2, (p1-p2).norm());
				m_constraints.push_back(c);
			}
		}
		break;

		case MESH_TYPE_PRIM:
		{
			// generate stretch constraints. assign a stretch constraint for each edge.
			EigenVector3 p1, p2;
			for (std::vector<Edge>::iterator e = m_mesh->m_edge_list.begin(); e != m_mesh->m_edge_list.end(); ++e)
			{
				p1 = m_mesh->m_current_positions.block_vector(e->m_v1);
				p2 = m_mesh->m_current_positions.block_vector(e->m_v2);
				SpringConstraint *c = new SpringConstraint(&m_stiffness_stretch, e->m_v1, e->m_v2, (p1 - p2).norm());
				m_constraints.push_back(c);
			}
		}
		break;

	}
	std::cout << m_constraints.size() << " LDALD " << std::endl;
}

void Simulation::calculateInertiaY()
{
	m_inertia_y = m_mesh->m_current_positions + m_mesh->m_current_velocities*0.9 * m_h;
}

void Simulation::calculateExternalForce()
{
    m_external_force.resize(m_mesh->m_system_dimension);
    m_external_force.setZero();

    // gravity
    for (unsigned int i = 0; i < m_mesh->m_vertices_number; ++i)
    {
        m_external_force[3*i+1] += -m_gravity_constant;
    }

    m_external_force = m_mesh->m_mass_matrix * m_external_force;
}

void Simulation::integrateClassicExplicitEuler()
{
	VectorX force;
	computeForces(m_mesh->m_current_positions, force);

	//vn+1 = vn + h*f(xn,vn)
	//xn+1 = xn + vn*h;
	m_mesh->m_previous_velocities = m_mesh->m_current_velocities;
	m_mesh->m_current_velocities = m_mesh->m_current_velocities + m_h*force;
	m_mesh->m_current_positions = m_mesh->m_current_positions + m_h* m_mesh->m_previous_velocities;

	//restaura posição
	
}

void Simulation::integrateExplicitEuler()
{
	VectorX force;
	computeForces(m_mesh->m_current_positions, force);

	VectorX pos_next;
	pos_next.resize(m_mesh->m_system_dimension);

	//two explicit euler different modes:
	if (m_integration_method == INTEGRATION_EXPLICIT_EULER_DISCRETE)
	{

		//update rule here is xn+1 = xn + m_h*(vn-1 + f(xn-1,vn-1))
		pos_next = m_mesh->m_current_positions +m_h*(m_mesh->m_previous_velocities + force*m_h);
	}
	else if (m_integration_method == INTEGRATION_EXPLICIT_EULER_SYMPLETIC)
	{
		//Sympletic Euler uses the new velocity when computing for the new position
		ScalarType mh_squared = m_h*m_h;

		//vn+1 = vn + h*f(xn,vn)
		//xn+1 = xn + vn+1*h;
		//update rule here is xn+1 = xn + m_h*(vn  + f(xn,vn))
		pos_next = m_mesh->m_current_positions + m_mesh->m_current_velocities*m_h + force*mh_squared;
	}

	updatePosAndVel(pos_next);
}

#pragma region implicit/explicit euler
void Simulation::computeForces(const VectorX& x, VectorX& force)
{
    VectorX gradient, dampened, resistance;

    gradient.resize(m_mesh->m_system_dimension);
    gradient.setZero();

	dampened.resize(m_mesh->m_system_dimension);
	dampened.setZero();

	resistance.resize(m_mesh->m_system_dimension);
	resistance.setZero();

    // springs
    for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
    {
        (*it)->EvaluateGradient(x, gradient);
    }

	// spring damping
	for (std::vector<Constraint*>::iterator it = m_constraints.begin(); it != m_constraints.end(); ++it)
	{
		(*it)->DampingGradient(x, m_mesh->m_current_velocities, dampened);
	}


	resistance = -m_mesh->m_current_velocities*(m_resistance_coefficient); //air resistance
	gradient += m_external_force;

    force = gradient + dampened + resistance;
}

#pragma endregion


#pragma region utilities
void Simulation::updatePosAndVel(const VectorX& new_pos)
{	
	m_mesh->m_previous_velocities = m_mesh->m_current_velocities;
    m_mesh->m_current_velocities = (new_pos - m_mesh->m_current_positions)/m_h;
    m_mesh->m_current_positions = new_pos;

	
	m_mesh->m_current_positions.block_vector(0) = GLM2Eigen(m_mesh->m_positions[0]);
	m_mesh->m_current_velocities.block_vector(0) = EigenVector3(0, 0, 0);

	m_mesh->m_current_positions.block_vector(59) = GLM2Eigen(m_mesh->m_positions[59]);
	m_mesh->m_current_velocities.block_vector(59) = EigenVector3(0, 0, 0);
}

#pragma endregion
