#ifndef _OPENGL_HEADERS_H_
#define _OPENGL_HEADERS_H_

// opengl headers
#include <GL/glew.h>
#include <GL/freeglut.h> 
#include <iostream>

#pragma comment(lib,"glew32.lib")
struct VBO
{
    VBO() //construtor
    {
        if(!glIsBuffer(m_vbo))
            glGenBuffers(1, &m_vbo); //gero um espaço no buffer para coordenadas xyz de cada vertice
        if(!glIsBuffer(m_cbo))
            glGenBuffers(1, &m_cbo); //gero um espaço no buffer para as cores rgb de cada vertice
        if(!glIsBuffer(m_nbo))
            glGenBuffers(1, &m_nbo); //gero um espaço no buffer para as normais rgb de cada vertice
        if(!glIsBuffer(m_tbo))
            glGenBuffers(1, &m_tbo); //gero um espaço no buffer para coordenadas de textura cada vertice
        if(!glIsBuffer(m_ibo))
            glGenBuffers(1, &m_ibo); //gero um espaço no buffer para os indices dos vertices que compoem os triangulos do modelo
    }

    virtual ~VBO() //desconstrutor, deleta as coisas do buffer
    {
        if(glIsBuffer(m_vbo))
            glDeleteBuffers(1, &m_vbo);
        if(glIsBuffer(m_cbo))
            glDeleteBuffers(1, &m_cbo);
        if(glIsBuffer(m_nbo))
            glDeleteBuffers(1, &m_nbo);
        if(glIsBuffer(m_tbo))
            glDeleteBuffers(1, &m_tbo);
        if(glIsBuffer(m_ibo))
            glDeleteBuffers(1, &m_ibo);
    }

    // vertex, color, normal, texture, triangle indeces
	// a funcao genBuffer associa as variaveis abaixo a buffers do opengl. O nome de cada um deles é um int.
    GLuint m_vbo, m_cbo, m_nbo, m_tbo, m_ibo;

    //0: modelview; 1: projection; 2: transformation; 3: enable_texture; 4: texture_sampler
	//essas são variaveis que a gente vai passar pro opengl la nos shaders. 
	//Similares aos buffers, os inteiros atrelam elas como nomes, a areas de memoria

    GLuint m_uniform_modelview;
    GLuint m_uniform_projection;
    GLuint m_uniform_transformation;
    GLuint m_uniform_enable_texture;
    GLuint m_uniform_texture_sampler;
	GLuint m_uniform_texture_sampler_2D;
	GLuint m_uniform_texture_sampler_tff;
};

#endif