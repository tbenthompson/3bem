#include "UnitTest++.h"
#include "mesh_3d.h"
#include <GL/glut.h>

std::vector<std::array<double, 3>> vertices = {
    {0.0,0.0,0.0},
    {0.0,0.0,1.0},
    {0.0,1.0,1.0},
    {1.0,1.0,0.0},
    {0.0,0.0,0.0},
    {0.0,1.0,0.0},
    {1.0,0.0,1.0},
    {0.0,0.0,0.0},
    {1.0,0.0,0.0},
    {1.0,1.0,0.0},
    {1.0,0.0,0.0},
    {0.0,0.0,0.0},
    {0.0,0.0,0.0},
    {0.0,1.0,1.0},
    {0.0,1.0,0.0},
    {1.0,0.0,1.0},
    {0.0,0.0,1.0},
    {0.0,0.0,0.0},
    {0.0,1.0,1.0},
    {0.0,0.0,1.0},
    {1.0,0.0,1.0},
    {1.0,1.0,1.0},
    {1.0,0.0,0.0},
    {1.0,1.0,0.0},
    {1.0,0.0,0.0},
    {1.0,1.0,1.0},
    {1.0,0.0,1.0},
    {1.0,1.0,1.0},
    {1.0,1.0,0.0},
    {0.0,1.0,0.0},
    {1.0,1.0,1.0},
    {0.0,1.0,0.0},
    {0.0,1.0,1.0},
    {1.0,1.0,1.0},
    {0.0,1.0,1.0},
    {1.0,0.0,1.0}
};

Mesh3D cube;
double rotate_y = 0; 
double rotate_x = 0;
void specialKeys( int key, int x, int y ) 
{
    if (key == GLUT_KEY_RIGHT)
        rotate_y += 5;
    if (key == GLUT_KEY_LEFT)
        rotate_y -= 5;
    if (key == GLUT_KEY_UP)
        rotate_x += 5;
    if (key == GLUT_KEY_DOWN)
        rotate_x -= 5;
    glutPostRedisplay();
}

void display()
{
    glClearColor( 0, 0, 0, 1 );
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    int w = glutGet( GLUT_WINDOW_WIDTH );
    int h = glutGet( GLUT_WINDOW_HEIGHT );
    gluPerspective( 60, w / h, 0.1, 100 );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt
        ( 
        3, 3, 3, 
        0, 0, 0,
        0, 0, 1
        );

    glRotatef( rotate_x, 1.0, 0.0, 0.0 );
    glRotatef( rotate_y, 0.0, 1.0, 0.0 );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    draw_mesh(cube);

    glutSwapBuffers();
}

int main(int argc, char **argv) {
    std::vector<std::array<int, 3>> faces;
    for (int i = 0; i < 12; i++) {
        faces.push_back({3 * i, 3 * i + 1, 3 * i + 2});
    }

    cube = {vertices, faces};
    std::cout << "Raw mesh has " << cube.vertices.size() << " vertices." << std::endl;
    cube = clean_mesh(cube);
    std::cout << "Cleaned mesh has " << cube.vertices.size() << " vertices." << std::endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(640, 480);
    glutCreateWindow("GLUT");
    glutDisplayFunc(display);
    glutSpecialFunc(specialKeys);
    glEnable(GL_DEPTH_TEST);
    // glutFullScreen();
    glutMainLoop();
}
