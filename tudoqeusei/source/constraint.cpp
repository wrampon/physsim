#include "constraint.h"

//----------Constraint Class----------//
Constraint::Constraint(ScalarType *stiffness) : 
    m_stiffness(stiffness)
{
}

Constraint::Constraint(const Constraint& other) : 
    m_stiffness(other.m_stiffness)
{
}

Constraint::~Constraint()
{
}

//----------AttachmentConstraint Class----------//
AttachmentConstraint::AttachmentConstraint(ScalarType *stiffness) : 
    Constraint(stiffness)
{
    m_selected = false;
}

AttachmentConstraint::AttachmentConstraint(ScalarType *stiffness, unsigned int p0, const EigenVector3& fixedpoint) : 
    Constraint(stiffness),
    m_p0(p0),
    m_fixd_point(fixedpoint)
{
    m_selected = false;
}

AttachmentConstraint::AttachmentConstraint(const AttachmentConstraint& other) : 
    Constraint(other),
    m_p0(other.m_p0),
    m_fixd_point(other.m_fixd_point),
    m_selected(other.m_selected)
{
    
}

AttachmentConstraint::~AttachmentConstraint()
{
}



// attachment spring gradient: k*(current_length)*current_direction
void AttachmentConstraint::EvaluateGradient(const VectorX& x, VectorX& gradient)
{
    EigenVector3 g_i = (*(m_stiffness))*(x.block_vector(m_p0) - m_fixd_point);
    gradient.block_vector(m_p0) += g_i;
}


void AttachmentConstraint::Draw(const VBO& vbos)
{
    m_attachment_constraint_body.move_to(Eigen2GLM(m_fixd_point));
    if (m_selected)
        m_attachment_constraint_body.change_color(glm::vec3(0.8, 0.8, 0.2));
    else
        m_attachment_constraint_body.change_color(glm::vec3(0.8, 0.2, 0.2));
        
    m_attachment_constraint_body.Draw(vbos);
}

//----------SpringConstraint Class----------//
SpringConstraint::SpringConstraint(ScalarType *stiffness) : 
    Constraint(stiffness)
{
}

SpringConstraint::SpringConstraint(ScalarType *stiffness, unsigned int p1, unsigned int p2, ScalarType length) : 
    Constraint(stiffness),
    m_p1(p1),
    m_p2(p2),
    m_rest_length(length)
{
}

SpringConstraint::SpringConstraint(const SpringConstraint& other) : 
    Constraint(other),
    m_p1(other.m_p1),
    m_p2(other.m_p2),
    m_rest_length(other.m_rest_length)
{
}

SpringConstraint::~SpringConstraint()
{
}


// sping gradient: k*(current_length-rest_length)*current_direction;
void SpringConstraint::EvaluateGradient(const VectorX& x, VectorX& gradient)
{
    EigenVector3 x_ij = x.block_vector(m_p2) - x.block_vector(m_p1);
    EigenVector3 g_ij = (*(m_stiffness))*(x_ij.norm()-m_rest_length)*x_ij.normalized()/2.0f;
    gradient.block_vector(m_p1) += (g_ij);
	gradient.block_vector(m_p2) -= (g_ij);
}

void SpringConstraint::DampingGradient(const VectorX& x, const VectorX& v ,VectorX& damping)
{
	EigenVector3 x_ij = x.block_vector(m_p2) - x.block_vector(m_p1);
	EigenVector3 v_ij = v.block_vector(m_p2) - v.block_vector(m_p1);
	EigenVector3 g_ij = x_ij.normalized()*(v_ij.dot(x_ij.normalized()));

	damping.block_vector(m_p1) += (g_ij);
	damping.block_vector(m_p2) -= (g_ij);
}