#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <FreeImage.h>
#include <vector>
#include <iostream>
#include <cmath>
// #include <Eigen/Dense>

#include "lagrangian.h"
#include "cloth.h"

#define mp 0.3f
#define STRECH_SPRING_L 5.0f
#define BEND_SPRING_L (STRECH_SPRING_L)
#define STRECH_SPRING_K 80.0f
#define BEND_SPRING_K 10.0f
#define MESH_SIZE 20
#define ACC_G glm::vec3(-2.0f,-2.0f,0.0f)
#define SPRING_DAMPING_COEFFICIENT 0.0f
#define TIME_STEP 0.01f

GLuint textureObj;

class System1 : public Lagrangian{
    public:
        int DOF = 3;
        glm::vec3 F{2.0f, 2.0f, 0.0f};
        float k = 80.0f;
        float k1 = 80.0f;
        float m = 0.3f;
        float l = 5.0f;
        float h = 0.01f;

        System1(glm::vec3 q, glm::vec3 dq, glm::vec3 p,std::vector<std::vector<particles>> sParticles): Lagrangian(q, dq, p) {

        }

        float L_evaluate(const glm::vec3& q, const glm::vec3& dq) {
            return 0.5*m*(dq.x*dq.x + dq.y*dq.y + dq.z*dq.z) - 0.5*k1*pow((sqrt((q.x*q.x + q.y*q.y + q.z*q.z)) - l) , 2) - (F.x*q.x + F.y*q.y + F.z*q.z);
        }     

        glm::vec3 dldq(const glm::vec3& q, const glm::vec3& dq) {
            glm::vec3 c;
            c.x = -k*q.x - F.x;
            c.y = -k*q.y - F.y;
            c.z = -k*q.z - F.z;
            return c;
        }

        glm::vec3 dldq_2(const glm::vec3& q, const glm::vec3& dq){
            glm::vec3 c;
            c.x = m*dq.x;
            c.y = m*dq.y;
            c.z = m*dq.z;
            return c;
        }

        std::vector<float> solveD1Ld(const std::vector<float>& q, const std::vector<float>& p, const std::vector<std::vector<particles>>& sParticles){
            std::vector<float> fq;
            for(int k=0; k < q.size(); k++){
                fq.push_back(0);
            }
            fq[0] = ((((p[0]/h) - (F[0]*m/2))+((k1/2))+ q[0]*(4.0f*m - k*h*h) + 2.0f*k*h*h*l)/(4.0f*m+k*h*h));
            return fq;
        }

        std::vector<float> solveD2Ld(const std::vector<float>& q, const std::vector<float>& fq, const std::vector<std::vector<particles>>& sParticles){
            std::vector<float> fp;
            for(int i=0; i < q.size(); i++){
                fp.push_back(0);
            }
            fp[0] = -(((2.0f*h*h*F[0]*m)+ k*h*h*(q[0]+fq[0]-2.0f*l) - 4.0f*m*(-q[0]+fq[0]))/(4.0f*h));
            return fp;
        }        

        void run(){
            
        }
};

clothsim::clothsim() {
    // particles
    for (int i = 0; i < 1; i++) {
        std::vector<particles> part;
        glm::vec3 origin(0.0f, 0.0f, 0.0f);
        for (int j = 0; j < MESH_SIZE; j++) {
            glm::vec2 texcoord(j / (MESH_SIZE - 1.0f), i / (MESH_SIZE - 1.0f));
            //glm::vec2 texcoord(0.0f, 0.0f);
            part.push_back(particles(mp, origin, texcoord));
            origin += glm::vec3(STRECH_SPRING_L, 0.0f, 0.0f);
        }
        vParticles.push_back(part);
    }

    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < MESH_SIZE - 1; j++)
            vStrechsprings.push_back(spring(&vParticles[i][j], &vParticles[i][j + 1], STRECH_SPRING_L, STRECH_SPRING_K));
    }
}


void clothsim::draw() {
    // Clear the color buffer and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Draw springs
    for (int i = 0; i < vStrechsprings.size(); i++) {
        vStrechsprings[i].draw();
    }

    // Draw particles
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < MESH_SIZE; j++) {
            vParticles[i][j].draw();
        }
    }
}


void clothsim::simulate(){
    //Clear the force accumulator for all particles
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < MESH_SIZE; j++)
            vParticles[i][j].removeforces();
    }

    lagrangeintegrator();
}


void clothsim::lagrangeintegrator(){
    float h =TIME_STEP;
    float l = 5.0f;
    for (int i = 0; i < 1; i++) {
        for (int j = 1; j < MESH_SIZE; j++) {
            particles &p = vParticles[i][j];

            glm::vec3 q = p.getposition();

            glm::vec3 dq = p.getvelocity();

            glm::vec3 pp = p.getmomentum();

            // System1 sy(q, dq, pp, vParticles);

            // Lagrangian* sys = &sy;

            if(j==1){
                float a1 = (p.getposition().x-vParticles[i][j-1].getposition().x);
                float b1 = (p.getposition().y-vParticles[i][j-1].getposition().y);
                float c1 = sqrt(a1*a1 + b1*b1);

                float a2 = -(p.getposition().x-vParticles[i][j+1].getposition().x);
                float b2 = -(p.getposition().y-vParticles[i][j+1].getposition().y);
                float c2 = sqrt(a2*a2 + b2*b2);

                float x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].getposition().x+vParticles[i][j+1].getposition().x+l*(a1/c1)-l*(a2/c2)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].getposition().x)+2.0f*(vParticles[i][j-1].getposition().x)-(vParticles[i][j+2].getposition().x)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].getposition().y+vParticles[i][j+1].getposition().y+l*(b1/c1)-l*(b2/c2)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].getposition().y)+2.0f*(vParticles[i][j-1].getposition().y)-(vParticles[i][j+2].getposition().y)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float z = 0.0f;

                p.setposition(glm::vec3(x, y, z));

                a1 = (p.getposition().x-vParticles[i][j-1].getposition().x);
                b1 = (p.getposition().y-vParticles[i][j-1].getposition().y);
                c1 = sqrt(a1*a1 + b1*b1);

                a2 = -(p.getposition().x-vParticles[i][j+1].getposition().x);
                b2 = -(p.getposition().y-vParticles[i][j+1].getposition().y);
                c2 = sqrt(a2*a2 + b2*b2);

                x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.getposition().x+q[0]-(vParticles[i][j-1].getposition().x+vParticles[i][j+1].getposition().x+l*(a1/c1)-l*(a2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.getposition().x+q[0])+vParticles[i][j+2].getposition().x-4.0f*(vParticles[i][j+1].getposition().x)+2.0f*(vParticles[i][j-1].getposition().x))-(mp/(h*h))*(p.getposition().x-q[0]));
                y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.getposition().y+q[1]-(vParticles[i][j-1].getposition().y+vParticles[i][j+1].getposition().y+l*(b1/c1)-l*(b2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.getposition().y+q[1])+vParticles[i][j+2].getposition().y-4.0f*(vParticles[i][j+1].getposition().y)+2.0f*(vParticles[i][j-1].getposition().y))-(mp/(h*h))*(p.getposition().y-q[1]));
                z = 0.0f;
                
                p.setmomentum(glm::vec3(x, y, z));
            }
            else if (j > 1 && j < MESH_SIZE - 2){
                float a1 = (p.getposition().x-vParticles[i][j-1].getposition().x);
                float b1 = (p.getposition().y-vParticles[i][j-1].getposition().y);
                float c1 = sqrt(a1*a1 + b1*b1);

                float a2 = -(p.getposition().x-vParticles[i][j+1].getposition().x);
                float b2 = -(p.getposition().y-vParticles[i][j+1].getposition().y);
                float c2 = sqrt(a2*a2 + b2*b2);

                float x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].getposition().x+vParticles[i][j+1].getposition().x+l*(a1/c1)-l*(a2/c2)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].getposition().x+vParticles[i][j-1].getposition().x)-(vParticles[i][j-2].getposition().x+vParticles[i][j+2].getposition().x)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].getposition().y+vParticles[i][j+1].getposition().y+l*(b1/c1)-l*(b2/c2)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].getposition().y+vParticles[i][j-1].getposition().y)-(vParticles[i][j-2].getposition().y+vParticles[i][j+2].getposition().y)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float z = 0.0f;

                p.setposition(glm::vec3(x, y, z));

                a1 = (p.getposition().x-vParticles[i][j-1].getposition().x);
                b1 = (p.getposition().y-vParticles[i][j-1].getposition().y);
                c1 = sqrt(a1*a1 + b1*b1);

                a2 = -(p.getposition().x-vParticles[i][j+1].getposition().x);
                b2 = -(p.getposition().y-vParticles[i][j+1].getposition().y);
                c2 = sqrt(a2*a2 + b2*b2); 

                x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.getposition().x+q[0]-(vParticles[i][j-1].getposition().x+vParticles[i][j+1].getposition().x+l*(a1/c1)-l*(a2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(3.0f*(p.getposition().x+q[0])+vParticles[i][j-2].getposition().x+vParticles[i][j+2].getposition().x-4.0f*(vParticles[i][j+1].getposition().x+vParticles[i][j-1].getposition().x))-(mp/(h*h))*(p.getposition().x-q[0]));
                y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.getposition().y+q[1]-(vParticles[i][j-1].getposition().y+vParticles[i][j+1].getposition().y+l*(b1/c1)-l*(b2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(3.0f*(p.getposition().y+q[1])+vParticles[i][j-2].getposition().y+vParticles[i][j+2].getposition().y-4.0f*(vParticles[i][j+1].getposition().y+vParticles[i][j-1].getposition().y))-(mp/(h*h))*(p.getposition().y-q[1]));
                z = 0.0f;

                p.setmomentum(glm::vec3(x, y, z));

            }
            else if(j==MESH_SIZE - 2){
                float a1 = (p.getposition().x-vParticles[i][j-1].getposition().x);
                float b1 = (p.getposition().y-vParticles[i][j-1].getposition().y);
                float c1 = sqrt(a1*a1 + b1*b1);

                float a2 = -(p.getposition().x-vParticles[i][j+1].getposition().x);
                float b2 = -(p.getposition().y-vParticles[i][j+1].getposition().y);
                float c2 = sqrt(a2*a2 + b2*b2);

                float x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].getposition().x+vParticles[i][j+1].getposition().x+l*(a1/c1)-l*(a2/c2)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j+1].getposition().x)+4.0f*(vParticles[i][j-1].getposition().x)-(vParticles[i][j-2].getposition().x)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].getposition().y+vParticles[i][j+1].getposition().y+l*(b1/c1)-l*(b2/c2)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j+1].getposition().y)+4.0f*(vParticles[i][j-1].getposition().y)-(vParticles[i][j-2].getposition().y)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float z = 0.0f;

                p.setposition(glm::vec3(x, y, z));

                a1 = (p.getposition().x-vParticles[i][j-1].getposition().x);
                b1 = (p.getposition().y-vParticles[i][j-1].getposition().y);
                c1 = sqrt(a1*a1 + b1*b1);

                a2 = -(p.getposition().x-vParticles[i][j+1].getposition().x);
                b2 = -(p.getposition().y-vParticles[i][j+1].getposition().y);
                c2 = sqrt(a2*a2 + b2*b2);

                x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.getposition().x+q[0]-(vParticles[i][j-1].getposition().x+vParticles[i][j+1].getposition().x+l*(a1/c1)-l*(a2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.getposition().x+q[0])+vParticles[i][j-2].getposition().x-2.0f*(vParticles[i][j+1].getposition().x)-4.0f*(vParticles[i][j-1].getposition().x))-(mp/(h*h))*(p.getposition().x-q[0]));
                y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.getposition().y+q[1]-(vParticles[i][j-1].getposition().y+vParticles[i][j+1].getposition().y+l*(b1/c1)-l*(b2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.getposition().y+q[1])+vParticles[i][j-2].getposition().y-2.0f*(vParticles[i][j+1].getposition().y)-4.0f*(vParticles[i][j-1].getposition().y))-(mp/(h*h))*(p.getposition().y-q[1]));
                z = 0.0f;

                p.setmomentum(glm::vec3(x, y, z));

            }
            else if (j==MESH_SIZE - 1){
                float a = (p.getposition().x-vParticles[i][j-1].getposition().x);
                float b = (p.getposition().y-vParticles[i][j-1].getposition().y);
                float c = sqrt(a*a + b*b);  
                float x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/4.0f)*(2.0f*vParticles[i][j-1].getposition().x+2.0f*l*(a/c)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/4.0f)-(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j-1].getposition().x)-(vParticles[i][j-2].getposition().x)))/((mp/(h*h)+(STRECH_SPRING_K/4.0f)+(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/4.0f)*(2.0f*vParticles[i][j-1].getposition().y+2.0f*l*(b/c)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/4.0f)-(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j-1].getposition().y)-(vParticles[i][j-2].getposition().y)))/((mp/(h*h)+(STRECH_SPRING_K/4.0f)+(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                float z = 0.0f;

                p.setposition(glm::vec3(x, y, z));

                a = (p.getposition().x-vParticles[i][j-1].getposition().x);
                b = (p.getposition().y-vParticles[i][j-1].getposition().y);
                c = sqrt(a*a + b*b);

                x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/4.0f)*(p.getposition().x+q[0]-2.0f*(vParticles[i][j-1].getposition().x+l*(a/c)))+(BEND_SPRING_K/(2.0f*l*l*l))*(0.5f*(p.getposition().x+q[0])+vParticles[i][j-2].getposition().x-2.0f*(vParticles[i][j-1].getposition().x))-(mp/(h*h))*(p.getposition().x-q[0]));
                y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/4.0f)*(p.getposition().y+q[1]-2.0f*(vParticles[i][j-1].getposition().y+l*(b/c)))+(BEND_SPRING_K/(2.0f*l*l*l))*(0.5f*(p.getposition().y+q[1])+vParticles[i][j-2].getposition().y-2.0f*(vParticles[i][j-1].getposition().y))-(mp/(h*h))*(p.getposition().y-q[1]));
                z = 0.0f;

                p.setmomentum(glm::vec3(x, y, z));

            }
        }
    }
}


clothsim cloth;

void DisplayCallback(void) {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    float center = (MESH_SIZE - 1) * STRECH_SPRING_L / 2.0f;
	gluLookAt(0.0f, 0.0f, center * 5, 0, 0, 0.0f, 0.0f, 1.0f, 0.0f); // Look at the center of the cloth


    // Draw the cloth simulation
    cloth.draw();

    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50.0f, 2.0f, 1.0f, center * 100);
	glutSwapBuffers();
}

// void SetTextures(const char *fileName) {
    
// }
// void SetLightSource(void) {
   
// }

// void SetMaterial(void) {
    
// }

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(900, 700);
    glutCreateWindow("Cloth Simulation");

    glutDisplayFunc(DisplayCallback);
    glutIdleFunc([]() { cloth.simulate(); glutPostRedisplay(); });

    glutMainLoop();
    return 0;
}
