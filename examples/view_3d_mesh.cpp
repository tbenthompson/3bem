#include "UnitTest++.h"
#include "mesh.h"
#include "mesh_gen.h"
#include "numerics.h"
#include <GL/glut.h>

void draw_mesh(Mesh& msh) {
    glBegin(GL_TRIANGLES);
    for (unsigned int i = 0; i < msh.faces.size(); i++) {
        for (int v = 0; v < 3; v++) {
            int vert = msh.faces[i][v];
            glVertex3f(msh.vertices[vert][0],
                       msh.vertices[vert][1],
                       msh.vertices[vert][2]);
        }
    }
    glEnd();
}

Mesh sphere;
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
    sphere = sphere_mesh({0.0, 1.0, 1.0}, 1.0);
    std::cout << "Raw mesh has " << sphere.vertices.size() << " vertices." << std::endl;
    sphere = refine_mesh(sphere, naturals(sphere.faces.size()));
    sphere = refine_mesh(sphere, naturals(sphere.faces.size()));
    sphere = refine_mesh(sphere, naturals(sphere.faces.size()));
    std::cout << "Refined mesh has " << sphere.vertices.size() << " vertices." << std::endl;
    sphere = clean_mesh(sphere);
    std::cout << "Cleaned mesh has " << sphere.vertices.size() << " vertices." << std::endl;

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
