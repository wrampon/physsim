#include <iostream>
#include "glsl_wrapper.h"
#include "stb_image.h"


RenderWrapper::RenderWrapper()
{
	
}
RenderWrapper::~RenderWrapper()
{
}

void RenderWrapper::InitShader(const char* vert_path, const char* frag_path)
{
    //Vamos a parte meio chatinha ^^

	//essas variaveis são inteiros, a funcao "glCreateShader" me retorna uma 
	//um inteiro que vai ser o "nome" daquela area de shader 
    m_vert_handle = glCreateShader(GL_VERTEX_SHADER);
    m_frag_handle = glCreateShader(GL_FRAGMENT_SHADER);

	//um shader eh um programa. Caso eu queira usar mais de um shader, preciso ter mais de um programa.
    m_shaderprog_handle = glCreateProgram();

    // carrego os arquivos contendo o texto do vertex e fragment shader
    const char* vert_source = textFileRead(vert_path);
    const char* frag_source = textFileRead(frag_path);

	//linko o texto com os "nomes" do vertex e fragment shader
    glShaderSource(m_vert_handle, 1, &vert_source, NULL);
    glShaderSource(m_frag_handle, 1, &frag_source, NULL);

	//compilo os shaders
    glCompileShader(m_vert_handle);
    glCompileShader(m_frag_handle);

    // verifico...ta tudo bem na compilação? 
    GLint compiled;
    glGetShaderiv(m_vert_handle, GL_COMPILE_STATUS, &compiled);
    if(!compiled)
        printShaderInfoLog(m_vert_handle);
    glGetShaderiv(m_frag_handle, GL_COMPILE_STATUS, &compiled);
    if(!compiled)
        printShaderInfoLog(m_frag_handle);

	//shaders nao sao debugaveis ): infelizmente

    //ligo os shaders lidos a um programa
    glAttachShader(m_shaderprog_handle, m_vert_handle);
    glAttachShader(m_shaderprog_handle, m_frag_handle);

	#if USE_GLSL_BEST_PRATICES == 0
		//má pratica informar no opengl os locations de cada coisa. 
		//as melhores praticas modernas dizem que tu fala o location dentro do shader direto
		//com a tag layout(location = n) em vez de fazer com o bind attrib location q nem aqui em baixo
		glBindAttribLocation(m_shaderprog_handle, VERTEX_COORD_ATTRIB, "v_position");
		glBindAttribLocation(m_shaderprog_handle, COLOR_ATRIB, "v_color");
		glBindAttribLocation(m_shaderprog_handle, NORMAL_ATTRIB, "v_normal");
		glBindAttribLocation(m_shaderprog_handle, TEXTURE_COORD_ATTRIB, "v_texcoord");

	#endif

	//linko o programa
    glLinkProgram(m_shaderprog_handle);
    GLint linked;
    //verifico se a linkagem foi efetiva
	glGetProgramiv(m_shaderprog_handle, GL_LINK_STATUS, &linked);
    if(!linked)
        printLinkInfoLog(m_shaderprog_handle);

    //Preencho algumas das variaveis que eu uso no shader...isso eu explico via audio ...é mais complicadinho

	m_vbo_handle.m_uniform_modelview = glGetUniformLocation(m_shaderprog_handle, "u_modelviewMatrix");
	m_vbo_handle.m_uniform_projection = glGetUniformLocation(m_shaderprog_handle, "u_projMatrix");
	m_vbo_handle.m_uniform_transformation = glGetUniformLocation(m_shaderprog_handle, "u_transformMatrix");
	m_vbo_handle.m_uniform_enable_texture = glGetUniformLocation(m_shaderprog_handle, "u_choose_tex");
	m_vbo_handle.m_uniform_texture_sampler = glGetUniformLocation(m_shaderprog_handle, "u_sampler1");

    // activate the shader program.
    glUseProgram(m_shaderprog_handle);
}
bool RenderWrapper::InitTexture(const char* tex_path)
{
    int x,y,n;
    unsigned char *data = stbi_load(tex_path, &x, &y, &n, 0);
	//essa biblioteca transforma uma png em um array de cores (tira os encodings e compactação)

    if (data == NULL) {
        fprintf(stderr, "Nao achei a textura.\n");
        return false;
    } else {
        GLint mode;
        if (n == 4) 
        {
            mode = GL_RGBA;
        }
        else if (n == 3) 
        {
            mode = GL_RGB;
        }
        else
        {
            fprintf(stderr, "Formato de textura invalido\n");
        }
        
		//aqui, to dizendo que a textura que eu vou usar agora é a que eu li
        glGenTextures(1, &m_texture);
		//é uma textura 2D (uma imagem)
        glBindTexture(GL_TEXTURE_2D, m_texture);

		//uma imagem, em modo rgb, com altura e largura x y, sem borda, com 32 bits, os dados tão em bytes, e estão nesse array data
		//basicamente, manda a textura para a placa de video
        glTexImage2D(GL_TEXTURE_2D, 0, mode, x, y, 0, mode, GL_UNSIGNED_BYTE, data);
        
		//tipo de interpolação que usa na textura (quando vista de perto e longe) pra deixar ela mais suave
		glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		//ativa o core de texturas 0
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(m_vbo_handle.m_uniform_enable_texture, 0);
		//linka a textura atual ao core 0
        glBindTexture(GL_TEXTURE_2D, m_texture);

		//a placa de video tem em torno de 8 cores de textura
		//caso tu queira usar mais de uma, tu pode só ativar sempre o mesmo, e mudar o bind_texture 
		//a utilidade de usar mais de um core é caso tu queira processar mais de uma textura no mesmo shader (explico isso em detalhes se tu quiser)

		//limpa o array
        free(data);
        return true;
    }
}
void RenderWrapper::CleanupShader()
{
	//desliga todo mundo xD
    glDetachShader(m_shaderprog_handle, m_vert_handle);
    glDetachShader(m_shaderprog_handle, m_frag_handle);
    glDeleteShader(m_vert_handle);
    glDeleteShader(m_frag_handle);
    glDeleteProgram(m_shaderprog_handle);
}

void RenderWrapper::SetCameraProjection(glm::mat4 projection)
{
	//bruxaria, ignora por enquanto
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(&projection[0][0]);
    
    ActivateShaderprog();
    glUniformMatrix4fv(m_vbo_handle.m_uniform_projection, 1, false, &projection[0][0]);
}
void RenderWrapper::SetCameraModelview(glm::mat4 modelview)
{
	//bruxaria, ignora por enquanto
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(&modelview[0][0]);
    
    ActivateShaderprog();
    glUniformMatrix4fv(m_vbo_handle.m_uniform_modelview, 1, false, &modelview[0][0]);
}

void RenderWrapper::ActivateShaderprog()
{
	//ativa o programa atual (o shader atual)
    GLint current_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &current_prog);
    if(current_prog != (GLint)m_shaderprog_handle)
        glUseProgram(m_shaderprog_handle);
}
void RenderWrapper::DeactivateShaderprog()
{
	//desativa o shader atual
    GLint current_prog;
    glGetIntegerv(GL_CURRENT_PROGRAM, &current_prog);
    if(current_prog == (GLint)m_shaderprog_handle)
        glUseProgram(0);
}

// funcoes auxiliares pra ler shader, e ver se compilou
char* RenderWrapper::textFileRead(const char* fileName) 
{
    char* text;

    if (fileName != NULL) {
        FILE *file = fopen(fileName, "rt");

        if (file != NULL) {
            fseek(file, 0, SEEK_END);
            int count = ftell(file);
            rewind(file);

            if (count > 0) {
                text = (char*)malloc(sizeof(char) * (count + 1));
                count = fread(text, sizeof(char), count, file);
                text[count] = '\0';    //cap off the string with a terminal symbol, fixed by Cory
            }
            fclose(file);
        }
    }
    return text;
}
void RenderWrapper::printLinkInfoLog(int prog) 
{
    int infoLogLen = 0;
    int charsWritten = 0;
    GLchar *infoLog;

    glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infoLogLen);

    // should additionally check for OpenGL errors here

    if (infoLogLen > 0)
    {
        infoLog = new GLchar[infoLogLen];
        // error check for fail to allocate memory omitted
        glGetProgramInfoLog(prog,infoLogLen, &charsWritten, infoLog);
        std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
        delete [] infoLog;
    }
}
void RenderWrapper::printShaderInfoLog(int shader)
{
    int infoLogLen = 0;
    int charsWritten = 0;
    GLchar *infoLog;

    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

    // should additionally check for OpenGL errors here

    if (infoLogLen > 0)
    {
        infoLog = new GLchar[infoLogLen];
        // error check for fail to allocate memory omitted
        glGetShaderInfoLog(shader,infoLogLen, &charsWritten, infoLog);
        std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
        delete [] infoLog;
    }

    // should additionally check for OpenGL errors here
}
