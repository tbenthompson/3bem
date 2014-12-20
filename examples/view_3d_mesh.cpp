#include "UnitTest++.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "numerics.h"
#include <GL/glut.h>

using namespace tbem;

void draw_mesh(Mesh<3> msh) {
    glBegin(GL_TRIANGLES);
    for(auto f: msh.facets) {
        for (auto v: f.vertices) {
            glVertex3f(v[0], v[1], v[2]);
        }
    }
    glEnd();
}

const Mesh<3> sphere = sphere_mesh({0.0, 1.0, 1.0}, 1.0).refine_repeatedly(3);
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
        0, 3, 0, 
        0, 0, 0,
        0, 0, 1
        );

    glRotatef( rotate_x, 1.0, 0.0, 0.0 );
    glRotatef( rotate_y, 0.0, 1.0, 0.0 );
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    draw_mesh(sphere);

    glutSwapBuffers();
}

#include "util.h"
int main(int argc, char **argv) {
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
