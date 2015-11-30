#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <vector>

#include "global_headers.h"
#include "mesh.h"
#include "constraint.h"


typedef enum
{
    INTEGRATION_EXPLICIT_EULER_COMMON,
	INTEGRATION_EXPLICIT_EULER_DISCRETE,
	INTEGRATION_EXPLICIT_EULER_SYMPLETIC
} IntegrationMethod;

class Simulation
{

public:
    Simulation();
	Simulation(IntegrationMethod im) :m_integration_method(im){}
    virtual ~Simulation();

    void Reset();
    void Update();
    void DrawConstraints(const VBO& vbos);


    inline void SetMesh(Mesh* mesh) {m_mesh = mesh;}
    
protected:

    // simulation constants
    ScalarType m_h; // time_step

    // simulation constants
    ScalarType m_gravity_constant;
    ScalarType m_stiffness_attachment;
    ScalarType m_stiffness_stretch;
    ScalarType m_stiffness_bending;
    ScalarType m_damping_coefficient;
	ScalarType m_resistance_coefficient;

    // integration method
    IntegrationMethod m_integration_method;

    // key simulation components: mesh and scene
    Mesh *m_mesh;
  
    // key simulation components: constraints
    std::vector<Constraint*> m_constraints;
    AttachmentConstraint* m_selected_attachment_constraint;

    // inertia term
    VectorX m_inertia_y;

    // external force (gravity, wind, etc...)
    VectorX m_external_force;

	bool g_gravity;


private:

    void clearConstraints(); // cleanup all constraints
    void setupConstraints(); // initialize constraints
    void calculateInertiaY(); // calculate the inertia term: y = current_pos + current_vel*h
    void calculateExternalForce(); // wind force is propotional to the area of triangles projected on the tangential plane

	void integrateClassicExplicitEuler();
	void integrateExplicitEuler();
	void PBD();

    // for explicit/implicit integration only
    void computeForces(const VectorX& x, VectorX& force);

    // evaluate gradient
    void evaluateGradient(const VectorX& x, VectorX& gradient);

    // utility functions
    void updatePosAndVel(const VectorX& new_pos); // current_vel = (next_pos-current_pos)/h; current_pos = next_pos; 

};

#endif