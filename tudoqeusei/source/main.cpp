#pragma warning( disable : 4244)
#include <iostream>
#include <string>

//----------Headers--------------//
#include "openGL_headers.h"
#include "math_headers.h"
#include "camera.h"
#include "primitive.h"
#include "glsl_wrapper.h"
#include "math.h"
#include "simulation.h"

//----------Global Parameters----------------//
int g_screen_width = DEFAULT_SCREEN_WIDTH;
int g_screen_height = DEFAULT_SCREEN_HEIGHT;
bool g_show_texture = true;
bool g_texture_load_succeed = false;
bool g_show_wireframe = false;	
bool first_time = true;
bool picker = false;
//----------Objects Parameters----------------//
int plane_side = 60; //lado do quadrilatero de tecido gerado
int g_max_fps = 60;
int g_time_step = 1000 / g_max_fps;
//----------Global Objects----------------//
Camera* g_camera;
std::vector<Mesh*> g_mesh_vector;
RenderWrapper* g_renderer;
Simulation* g_simulation;

//----------Mouse Control--------------------//
glm::vec2 g_mouse_old;
glm::vec2 dx_mouse;
glm::vec2 g_mouse_press;
int g_mouse_wheel_pos;	
int it = 0;
unsigned char g_button_mask = 0x00;
float mouseZ;

//----------glut function handlers-----------//
void resize(int, int);
void timeout(int);
void display(void);
void mouse_motion(int, int);
void drawGrid();
void mouse_click(int button, int state, int x, int y);
void CallBackKeyboardFunc(unsigned char, int, int);

//----------other utility functions----------//
void init(void);
glm::vec3 gluUnprojectUso();
//----------misc variables------------------//

int main(int argc, char ** argv)
{
    //A funcao main do openGL eh igual a main comum de c, chamada inicialmente
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA|GLUT_MULTISAMPLE);

	glutInitContextVersion (3 , 1);
	glutInitContextProfile (GLUT_CORE_PROFILE );

    glutCreateWindow("Licao infinita");

	//inicia a janelinha
    glutInitWindowSize(g_screen_width, g_screen_height);

    //inicializacao do openGL
    init(); //minha funcao init
    glutReshapeWindow(g_screen_width, g_screen_height);

    //gruda as funcoes basicas aos eventos:
		
	//sempre que for feito um refresh de tela pelo monitor, chamarei essa funcao
    glutDisplayFunc(display);

	//sempre que janela for redimensionada, chamarei essa funcao
    glutReshapeFunc(resize);

	//sempre que tiver tempo livre do processador, chamo essa funcao
    glutIdleFunc(display);

	//se uma tecla for apertada, chamo essa funcao
	glutKeyboardFunc(CallBackKeyboardFunc);

	//se o mouse for movimentado dentro da janela do openGL chamo essa funcao
	glutMouseFunc(mouse_click);
    glutMotionFunc(mouse_motion);

    glutMainLoop();
    return 0;
}


void init(void)
{
	//inicializando shaders
    fprintf(stdout, "Inicializando Shaders...\n");

	//seta a glew para usar o opengl novo
	glewExperimental = GL_TRUE;

	//glew é o cara que faz o papel de lidar com shaders
	if(glewInit()!= GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
        exit(EXIT_FAILURE);
    }

    printf ("Vendor: %s\n", glGetString (GL_VENDOR));
    printf ("Renderer: %s\n", glGetString (GL_RENDERER));
    printf ("Version: %s\n", glGetString (GL_VERSION));
    printf ("GLSL: %s\n", glGetString (GL_SHADING_LANGUAGE_VERSION));

	//criando uma nova classe que cuida do render e de ligar os shaders
    g_renderer = new RenderWrapper();

	//carregando os shaders padrão (entra pra ver detalhes)
    g_renderer->InitShader(DEFAULT_VERT_SHADER_FILE, DEFAULT_FRAG_SHADER_FILE);

	//carregando uma textura padrão (entra pra ver detalhes)
    g_texture_load_succeed = g_renderer->InitTexture(DEFAULT_TEXTURE_FILE);

	//criando um objeto camera (mistica)
	g_camera = new Camera();

	//redimensionando a camera pra usar o tamanho de janela padrao
	g_camera->Reset(g_screen_width, g_screen_height);

	//habilitando teste de profundidade
    glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glEnable(GL_MULTISAMPLE);

	//criando um objeto do tipo plano com lado = 100 (100 vertices de amostragem)
	Mesh* g_mesh;

	g_mesh = new Plane(plane_side);
	g_mesh_vector.push_back(g_mesh);
	//inicializando o plano
	g_mesh_vector[0]->Reset();

	g_mesh = new ObjMesh(DEFAULT_OBJ_FILE);
	g_mesh_vector.push_back(g_mesh);
	g_mesh_vector[1]->Reset();

	//cria simulação
	g_simulation = new Simulation(INTEGRATION_EXPLICIT_EULER_SYMPLETIC);
	g_simulation->SetMesh(g_mesh_vector[0]);
	g_simulation->Reset();
}
void display() 
{
	int initTime = glutGet(GLUT_ELAPSED_TIME);
	
	if ((g_button_mask & 0x08))
	{
		if (!picker){
			glm::vec3 mouse = gluUnprojectUso();
			g_mesh_vector[0]->grabMesh(mouse);
			picker = true;
		}
		g_mesh_vector[0]->calculatePointDisplacement(*g_camera, g_mouse_press, g_mouse_old, mouseZ);
		g_mouse_press = g_mouse_old;
	}
	else
	{
		picker = false;
	}
	//Limpa buffer de cores, e de teste de profundidade
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//Fundo branco
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

	//Miro a camera no universo (mais caixa preta, splico uma hora dessas)
	g_renderer->SetCameraModelview(g_camera->GetViewMatrix());
	g_renderer->SetCameraProjection(g_camera->GetProjectionMatrix());

	//um gerador de terreno tosco... altera os vertices pra mudar topologia...entrar pra ver
	//if(it<200)
	//g_mesh_vector[0]->displacePoint();

	//quantos steps da simulação eu faço antes de renderizar ?
	for (int i = 0; i < SUBSTEPS; i++)
	{
		if (picker)g_mesh_vector[0]->displacePoint();

		g_simulation->Update();
	}

	g_renderer->ActivateShaderprog(); //ativa shader atual
		//caso queira entender melhor, entra na funcao draw da Mesh
		g_mesh_vector[0]->Draw(g_renderer->getVBO(),g_show_wireframe,g_show_texture & g_texture_load_succeed); //desenha plano usando a gpu
		//g_mesh_vector[1]->Draw(g_renderer->getVBO(),g_show_wireframe,false & g_texture_load_succeed); //desenha coelho usando a gpu
	g_renderer->DeactivateShaderprog(); //desativa shader atual

	//g_camera->DrawAxis(); //desenha eixos da camera sem usar shader

	//Troco o buffer de desenho atual pelo outro buffer (tem dois buffers que se alternam)
    glutSwapBuffers();

	it++;
	//Sleep(maximun(1000.0f / 30.f - (glutGet(GLUT_ELAPSED_TIME) - initTime), 0));
}
void resize(int width, int height)
{    
	//mais um pouquinho de magica nesse ratio, explico em um momento futuro
	g_screen_width = width;
    g_screen_height = height;
    glViewport(0, 0, width, height);
	g_camera->ResizeWindow(width, height);
    glutPostRedisplay();
}
void keyboard ( unsigned char key, int x, int y ) {
	switch ( key ) {
	case 27:
		exit (0);
		break;
	default:
		break;
	}
}


void CallBackKeyboardFunc(unsigned char key, int x, int y)
{
	switch (key) {
	case 'w':
		g_mesh_vector[0]->m_current_positions.block_vector(0) += EigenVector3(0, 0.1, 0);
	break;

	case 'f':
		g_camera->targetPosition(glm::vec3(5, -3, 2), glm::vec3(0, 0, 0));
		break;
		
	}
}

glm::vec3 gluUnprojectUso()
{
	int viewport[4];
	double modelview[16];
	double projection[16];
	double posX, posY, posZ;	
	
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

	float wy = (float)viewport[3] - (float)g_mouse_old.y;

	float winX = g_mouse_old.x;
	float winY = (float)viewport[3] - (float)g_mouse_old.y;
	float winZ;
	glReadPixels((int)g_mouse_old.x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

	return glm::vec3(posX, posY, posZ);
}

void mouse_click(int button, int state, int x, int y)
{
	switch(state)
	{
		case GLUT_DOWN:
			if (glutGetModifiers() != GLUT_ACTIVE_CTRL)
			{
				// left: 0. right: 2. middle: 1.
				g_button_mask |= 0x01 << button;
				g_mouse_old.x = x;
				g_mouse_old.y = y;
				g_mouse_press.x = x;
				g_mouse_press.y = y;
			}
			else if (glutGetModifiers() == GLUT_ACTIVE_CTRL)
			{
				// ctrl: 3
				g_button_mask |= 0x01 << 3;
				g_mouse_old.x = x;
				g_mouse_old.y = y;
				g_mouse_press.x = x;
				g_mouse_press.y = y;
			}
			break;
		case GLUT_UP:
			if (glutGetModifiers() == GLUT_ACTIVE_CTRL)
			{// special case for ctrl
				button = 3;
			}
			unsigned char mask_not = ~g_button_mask;
			mask_not |= 0x01 << button;
			g_button_mask = ~mask_not;
			break;
	}

}

void mouse_motion(int x, int y)
{
	dx_mouse.x = (float)(x - g_mouse_old.x);
	dx_mouse.y = (float)(y - g_mouse_old.y);

	if (g_button_mask & 0x01)//botao esquerdo do mouse
    {
		g_camera->MouseChangeHeadPitch(0.2f, dx_mouse.x, dx_mouse.y); //rotaciona camera em volta do centro
    } 
    else if (g_button_mask & 0x02) //botao scroll apertado
    {
		g_camera->MouseChangeLookat(0.001f, dx_mouse.x, dx_mouse.y); //muda o centro de onde a camera ta olhando
    }
    else if (g_button_mask & 0x04) //botao direito do mouse
    {
		g_camera->MouseChangeDistance(0.05f, dx_mouse.x, dx_mouse.y);//muda zoom da camera
    }

    g_mouse_old.x = x;
    g_mouse_old.y = y;
   
}