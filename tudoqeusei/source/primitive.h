#ifndef _PRIMITIVE_H_
#define _PRIMITIVE_H_

#include "math_headers.h"
#include "opengl_headers.h"
#include "mesh.h"
   
enum PrimitiveType {PLANE, SPHERE, CUBE};

class Primitive : public Mesh
{
public:
    Primitive(const PrimitiveType& t) : Mesh(MESH_TYPE_PRIM), m_type(t), m_pos(glm::vec3(0,0,0)), m_vel(glm::vec3(0,0,0)),  m_has_vel(false), m_has_gravity(false), m_previous_pos(glm::vec3(0,0,0)) {};
   
    PrimitiveType type() const
    {
        return m_type;
    }   
	virtual void Draw2();
	virtual void generateTriangleList();
	virtual void generateContours(){std::cout<<"Warning: base class Primitive initialization reached"<<std::endl;}
	virtual bool Init() {std::cout<<"Warning: base class Primitive initialization reached"<<std::endl; return false;}
public:
    glm::vec3 m_pos;
    glm::vec3 m_previous_pos;
    glm::vec3 m_vel;
    bool m_has_vel;
    bool m_has_gravity;

protected:
    PrimitiveType m_type;

};

class Plane : public Primitive
{
public:
	Plane(): Primitive(PLANE), side(10){};
	Plane(int side): Primitive(PLANE), side(side){};
	virtual ~Plane(){};
	virtual bool Init();
	virtual void generateTriangleList();
	int side;
	void createMountains();
};

class Sphere : public Primitive
{
public:
	Sphere() : Primitive(SPHERE), m_radius(1){};
    Sphere(float radius) : Primitive(SPHERE), m_radius(radius) {};
    virtual ~Sphere() {};
	virtual bool Init();
	float m_radius;
	void generateEdgeList();
};

class Cube : public Primitive
{
public:
	//Construtores
	Cube(int side): Primitive(CUBE), side(side){};

    virtual ~Cube() {};
	virtual bool Init();
	void generateEdgeList();
	virtual void generateTriangleList();

protected:
	int side;
};


#endif