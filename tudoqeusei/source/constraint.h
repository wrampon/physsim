#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

#include <vector>
#include <iostream>

#include "global_headers.h"
#include "math_headers.h"
#include "opengl_headers.h"
#include "primitive.h"

class Constraint
{
public:
    Constraint(ScalarType *stiffness);
    Constraint(const Constraint& other);
    virtual ~Constraint();

    virtual ScalarType  EvaluatePotentialEnergy(const VectorX& x) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl; return 0;}
    virtual void  EvaluateGradient(const VectorX& x, VectorX& gradient) {std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl;}
	virtual void  DampingGradient(const VectorX& x, const VectorX& v, VectorX& gradient) { std::cout << "Warning: reach <Constraint> base class virtual function." << std::endl; }

    inline const ScalarType& Stiffness() {return (*m_stiffness);}

protected:
    ScalarType *m_stiffness;

// for visualization and selection
public:
    virtual void Draw(const VBO& vbos) { /*do nothing*/ }
    //virtual ScalarType RayConstraintIntersection() {return false;}
};

class AttachmentConstraint : public Constraint
{
public:
    AttachmentConstraint(ScalarType *stiffness);
    AttachmentConstraint(ScalarType *stiffness, unsigned int p0, const EigenVector3& fixedpoint);
    AttachmentConstraint(const AttachmentConstraint& other);
    virtual ~AttachmentConstraint();

    virtual void  EvaluateGradient(const VectorX& x, VectorX& gradient);
	virtual void  DampingGradient(const VectorX& x, const VectorX& v, VectorX& gradient){}

protected:
    unsigned int m_p0;
    EigenVector3 m_fixd_point;

// for visualization and selection
public:
    virtual void Draw(const VBO& vbos);
    inline void Select() {m_selected = true;}
    inline void UnSelect() {m_selected = false;}
    inline EigenVector3 GetFixedPoint() {return m_fixd_point;}
    inline void SetFixedPoint(const EigenVector3& target) {m_fixd_point = target;}
    inline unsigned int GetConstrainedVertexIndex() {return m_p0;}

private: 
    bool m_selected;
    Sphere m_attachment_constraint_body;
};

class SpringConstraint : public Constraint
{
public:
    SpringConstraint(ScalarType *stiffness);
    SpringConstraint(ScalarType *stiffness, unsigned int p1, unsigned int p2, ScalarType length);
    SpringConstraint(const SpringConstraint& other);
    virtual ~SpringConstraint();

    virtual void  EvaluateGradient(const VectorX& x, VectorX& gradient);
	virtual void  DampingGradient(const VectorX& x, const VectorX& v,VectorX& gradient);

protected:
    unsigned int m_p1, m_p2;
    // rest length
    ScalarType m_rest_length;
};

#endif