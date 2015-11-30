#ifndef _COMMON_HEADERS_H_
#define _COMMON_HEADERS_H_

// glsl
#define USE_GLSL_BEST_PRATICES 1

// precision
#define HIGH_PRECISION

// single or double presicion
#ifdef HIGH_PRECISION
    typedef double ScalarType;
    #define TW_TYPE_SCALAR_TYPE TW_TYPE_DOUBLE
#else
    typedef float ScalarType;
    #define TW_TYPE_SCALAR_TYPE TW_TYPE_FLOAT
#endif

// small number and large number
#ifdef HIGH_PRECISION
    #define EPSILON 1e-15
#else
    #define EPSILON 1e-6
#endif

#define LARGER_EPSILON 1e-6
#define EPSILON_GAMBI 0.000001

// default values
#define DEFAULT_SCREEN_WIDTH 1024
#define DEFAULT_SCREEN_HEIGHT 768
#define PI 3.14159265359
#define SUBSTEPS 10

// selection radius
#define DEFAULT_SELECTION_RADIUS 1.0

// file localtions
#define DEFAULT_VERT_SHADER_FILE "./shaders/vert.glsl"
#define DEFAULT_FRAG_SHADER_FILE "./shaders/frag.glsl"
#define DEFAULT_TEXTURE_FILE "./textures/checkerboard.png"
#define DEFAULT_OBJ_FILE "./obj/bunny.obj"

#endif