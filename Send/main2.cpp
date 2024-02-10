#include <FreeImage.h>
#include <iostream>
#include <cmath>

#include "time_integrator.h"
// #include "rope.h"
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>

#define mp 0.3f
#define STRECH_SPRING_L 5.0f
#define KS 80.0f
#define ROPE_SIZE 20

#define FORCE_G std::vector<float>{-0.0f,-1.8f,0.0f}
#define H 0.05f

class System : public Integrator{
private:
    std::vector<Eigen::VectorXf>  initialise_q(){
        // std::cout << "Hello q\n";
        std::vector<Eigen::VectorXf> a;
        for(int i = 0; i<ROPE_SIZE; i++){
            Eigen::VectorXf t(5);
            t << i*STRECH_SPRING_L, 0.0, 0.0, 0.0, 0.0;
            a.push_back(t);

            // for(int j = 0; j<5; j++){
            //     std::cout << a[i][j] << " ";
            // }
            // std::cout << std::endl;
        }

        // std::cout << "End q\n";

        return a;
    }

    std::vector<Eigen::VectorXf>  initialise_p(){
        std::vector<Eigen::VectorXf> a;
        // std::cout << "Hello p\n";
        for(int i = 0; i<ROPE_SIZE; i++){
            Eigen::VectorXf t(5);
            t << 0.0, 0.0, 0.0, 0.0, 0.0;
            a.push_back(t);

            // for(int j = 0; j<5; j++){
            //     std::cout << a[i][j] << " ";
            // }
            // std::cout << std::endl;
        }

        // std::cout << "End p\n";

        return a;
    }

public:
    // Every particle out there will have position, amplitude and phase

    System() : Integrator(initialise_q(), initialise_p()){

    }


    void update(){
        std::vector<Eigen::VectorXf> old_q = q;
        std::vector<Eigen::VectorXf> old_p = p;

        // for(int i = 1; i< 2; i++){
        //     std::cout << "Value of point " << i << " = { " << old_q[i][0] << " " << old_q[i][1] << " " << old_q[i][2] << " " << old_q[i][3] << " " << old_q[i][4] << " }\n";
        // }

        for(int i = 1; i< ROPE_SIZE; i++){

            if(i != ROPE_SIZE - 1){

                float dist_A = std::sqrt(std::pow(old_q[i][0] - old_q[i-1][0], 2) + std::pow(old_q[i][1] - old_q[i-1][1], 2) + std::pow(old_q[i][2] - old_q[i-1][2], 2));
                float dist_B = std::sqrt(std::pow(old_q[i][0] - old_q[i+1][0], 2) + std::pow(old_q[i][1] - old_q[i+1][1], 2) + std::pow(old_q[i][2] - old_q[i+1][2], 2));
                float delphi_A = std::pow(old_q[i][4] - old_q[i-1][4], 2);
                float delphi_B = std::pow(old_q[i][4] - old_q[i+1][4], 2); 
                float delamp_A = std::pow(old_q[i][3] - old_q[i-1][3], 2);
                float delamp_B = std::pow(old_q[i][3] - old_q[i+1][3], 2); 
                float camp_A = old_q[i-1][3]*old_q[i-1][3] + old_q[i-1][3]*old_q[i][3] + old_q[i][3]*old_q[i][3];
                float camp_B = old_q[i+1][3]*old_q[i+1][3] + old_q[i+1][3]*old_q[i][3] + old_q[i][3]*old_q[i][3];

                float const_A = -KS/2*(delphi_A*camp_A/6 + delamp_A/2)/std::pow(dist_A, 3);
                float const_B = -KS/2*(delphi_B*camp_B/6 + delamp_B/2)/std::pow(dist_B, 3);

                // std::cout << dist_A << " " << delamp_A << " " << delphi_A << " " << camp_A << " " << const_A << "\n";
                // std::cout << delamp_B << " " << delphi_B << " " << camp_B << " " << const_B << "\n";

                q[i][0] = ( -FORCE_G[0] + old_p[i][0]/H + old_q[i][0]*mp/(H*H) - KS*(dist_A - 1.0)*(old_p[i][0] - old_p[i-1][0])/dist_A -KS*(dist_B - 1.0)*(old_p[i][0] - old_p[i+1][0])/dist_B + const_A*(old_p[i][0] - old_p[i-1][0]) + const_B*(old_p[i][0] - old_p[i+1][0])) / ( mp / (H*H) );

                q[i][1] = ( -FORCE_G[1] + old_p[i][1]/H + old_q[i][1]*mp/(H*H) - KS*(dist_A - 1.0)*(old_p[i][1] - old_p[i-1][1])/dist_A -KS*(dist_B - 1.0)*(old_p[i][1] - old_p[i+1][1])/dist_B + const_A*(old_p[i][1] - old_p[i-1][1]) + const_B*(old_p[i][1] - old_p[i+1][1])) / ( mp / (H*H) );

                q[i][2] = ( -FORCE_G[2] + old_p[i][2]/H + old_q[i][2]*mp/(H*H) - KS*(dist_A - 1.0)*(old_p[i][2] - old_p[i-1][2])/dist_A -KS*(dist_B - 1.0)*(old_p[i][2] - old_p[i+1][2])/dist_B + const_A*(old_p[i][2] - old_p[i-1][2]) + const_B*(old_p[i][2] - old_p[i+1][2])) / ( mp / (H*H) );

                p[i][0] = H*(mp/(H*H)*(q[i][0] - old_q[i][0]) + FORCE_G[0] + KS*(dist_A - 1.0)*(old_p[i][0] - old_p[i-1][0])/dist_A  + KS*(dist_B - 1.0)*(old_p[i][0] - old_p[i+1][0])/dist_B - const_A*(old_p[i][0] - old_p[i-1][0]) - const_B*(old_p[i][0] - old_p[i+1][0]));

                p[i][1] = H*(mp/(H*H)*(q[i][1] - old_q[i][1]) + FORCE_G[1] + KS*(dist_A - 1.0)*(old_p[i][1] - old_p[i-1][1])/dist_A  + KS*(dist_B - 1.0)*(old_p[i][1] - old_p[i+1][1])/dist_B - const_A*(old_p[i][1] - old_p[i-1][1]) - const_B*(old_p[i][1] - old_p[i+1][1]));
            
                p[i][2] = H*(mp/(H*H)*(q[i][2] - old_q[i][2]) + FORCE_G[2] + KS*(dist_A - 1.0)*(old_p[i][2] - old_p[i-1][2])/dist_A  + KS*(dist_B - 1.0)*(old_p[i][2] - old_p[i+1][2])/dist_B - const_A*(old_p[i][2] - old_p[i-1][2]) - const_B*(old_p[i][2] - old_p[i+1][2]));
            }
            else{

                float dist_A = std::sqrt(std::pow(old_q[i][0] - old_q[i-1][0], 2) + std::pow(old_q[i][1] - old_q[i-1][1], 2) + std::pow(old_q[i][2] - old_q[i-1][2], 2));
                float delphi_A = std::pow(old_q[i][4] - old_q[i-1][4], 2);
                float delamp_A = std::pow(old_q[i][3] - old_q[i-1][3], 2);
                float camp_A = old_q[i-1][3]*old_q[i-1][3] + old_q[i-1][3]*old_q[i][3] + old_q[i][3]*old_q[i][3];


                float const_A = -KS/2*(delphi_A*camp_A/6 + delamp_A/2)/std::pow(dist_A, 3);

                q[i][0] = ( -FORCE_G[0] + old_p[i][0]/H + old_q[i][0]*mp/(2*H*H) - KS*(dist_A - 1.0)*(old_p[i][0] - old_p[i-1][0])/dist_A  + const_A*(old_p[i][0] - old_p[i-1][0]) ) / ( mp / (2*H*H) );

                q[i][1] = ( -FORCE_G[1] + old_p[i][1]/H + old_q[i][1]*mp/(2*H*H) - KS*(dist_A - 1.0)*(old_p[i][1] - old_p[i-1][1])/dist_A  + const_A*(old_p[i][1] - old_p[i-1][1]) ) / ( mp / (2*H*H) );

                q[i][2] = ( -FORCE_G[2] + old_p[i][2]/H + old_q[i][2]*mp/(2*H*H) - KS*(dist_A - 1.0)*(old_p[i][2] - old_p[i-1][2])/dist_A  + const_A*(old_p[i][2] - old_p[i-1][2]) ) / ( mp / (2*H*H) );

                p[i][0] = H*(mp/(2*H*H)*(q[i][0] - old_q[i][0]) + FORCE_G[0] + KS*(dist_A - 1.0)*(old_p[i][0] - old_p[i-1][0])/dist_A   - const_A*(old_p[i][0] - old_p[i-1][0]) );

                p[i][1] = H*(mp/(2*H*H)*(q[i][1] - old_q[i][1]) + FORCE_G[1] + KS*(dist_A - 1.0)*(old_p[i][1] - old_p[i-1][1])/dist_A   - const_A*(old_p[i][1] - old_p[i-1][1]) );
            
                p[i][2] = H*(mp/(2*H*H)*(q[i][2] - old_q[i][2]) + FORCE_G[2] + KS*(dist_A - 1.0)*(old_p[i][2] - old_p[i-1][2])/dist_A  - const_A*(old_p[i][2] - old_p[i-1][2]) );

            }

        }
        

    }


    void run(){
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        for(int i = 0; i< ROPE_SIZE; i++){
            glColor3f(1.0f, 1.0f, 1.0f);
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glVertex3fv(glm::value_ptr(glm::vec3{q[i][0], q[i][1], q[i][2]}));
            glEnd();
        }

        for (int i = 0; i < ROPE_SIZE - 1; i++) {
            glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_LINES);
            glVertex3fv(glm::value_ptr(glm::vec3{q[i][0], q[i][1], q[i][2]}));
            glVertex3fv(glm::value_ptr(glm::vec3{q[i+1][0], q[i+1][1], q[i+1][2]}));
            glEnd();    
        }
        

    }

    void simulate(){
        update();
    }



};

System rope;



void DisplayCallback(void) {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    float center = (ROPE_SIZE - 1) * STRECH_SPRING_L / 2.0f;
	gluLookAt(0.0f, 0.0f, center * 5, 0, 0, 0.0f, 0.0f, 1.0f, 0.0f); // Look at the center of the cloth


    // Draw the cloth simulation
    rope.run();

    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(50.0f, 2.0f, 1.0f, center * 100);
	glutSwapBuffers();
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(900, 700);
    glutCreateWindow("1D Simulation : Version 2");

    glutDisplayFunc(DisplayCallback);
    glutIdleFunc([]() { rope.simulate(); glutPostRedisplay(); });

    glutMainLoop();
    return 0;
}