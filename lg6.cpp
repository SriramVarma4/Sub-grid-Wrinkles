#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <FreeImage.h>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

#define mp 0.3f
#define STRECH_SPRING_L 5.0f
#define BEND_SPRING_L (STRECH_SPRING_L)
#define STRECH_SPRING_K 80.0f
#define BEND_SPRING_K 10.0f
#define MESH_SIZE 10
#define ACC_G glm::vec3(-2.0f,-2.0f,0.0f)
#define SPRING_DAMPING_COEFFICIENT 0.0f
#define TIME_STEP 0.01f

GLuint textureObj;



class particles {
public:
    particles(float, glm::vec3, glm::vec2);

    void draw(void);
    void removeforces(void);

    const float m_mass;
    glm::vec3 m_pos;
    glm::vec3 m_vel;
    glm::vec3 m_mom;
    glm::vec3 m_acc;
    glm::vec3 m_totalforce;

    glm::vec2 m_texcoord;
    glm::vec3 m_normal;
};

particles::particles(float mass, glm::vec3 position, glm::vec2 texcoord) : m_mass(mass) {
    m_pos = position;
    m_vel = glm::vec3(0.0f, 0.0f, 0.0f);
    m_mom = glm::vec3(0.0f, 0.0f, 0.0f);
    m_acc = glm::vec3(0.0f, 0.0f, 0.0f);
    m_totalforce = glm::vec3(0.0f, 0.0f, 0.0f);

    m_texcoord = texcoord;
    m_normal = glm::vec3(0.0f, 0.0f, 0.0f);
}

void particles::draw(void) {
    glColor3f(1.0f, 1.0f, 1.0f);
    glPointSize(3.0f);
    glBegin(GL_POINTS);
    glVertex3fv(glm::value_ptr(m_pos));
    glEnd();
}

void particles::removeforces(void) {
    m_totalforce = glm::vec3(0.0f, 0.0f, 0.0f);
}


class spring {
public:
    spring(particles *, particles *, float, float);

    void draw(void);
    float deflength(void);
    glm::vec3 springforce(void);

    particles *m_particle1;
    particles *m_particle2;
    const float m_L;
    const float m_K;
};

spring::spring(particles *p1, particles *p2, float L, float K) : m_L(L), m_K(K) {
    m_particle1 = p1;
    m_particle2 = p2;
}

void spring::draw(void) {
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    glVertex3fv(glm::value_ptr(m_particle1->m_pos));
    glVertex3fv(glm::value_ptr(m_particle2->m_pos));
    glEnd();
}

float spring::deflength(void) {
    glm::vec3 defor = m_particle1->m_pos - m_particle2->m_pos;
    return glm::length(defor);
}

glm::vec3 spring::springforce(void) {
    glm::vec3 link = m_particle2->m_pos - m_particle1->m_pos;
    float slength = glm::length(link);
    glm::vec3 force = (link - (link / slength) * m_L) * m_K;
    return force;
}


class L{
    vector<float> q_k;
    vector<float> dq_k;
    vector<float> p_k;

    public:
        L(vector<float> q, vector<float> dq, vector<float> p){
            q_k = q;
            dq_k = dq;
            p_k = p;
        }

        virtual float L_evaluate(vector<float> q, vector<float> dq){
            return 0.0;
        }

        virtual vector<float> dldq(vector<float> q, vector<float> dq){
            vector<float> c;
            return c;
        }

        virtual vector<float> dldq_2(vector<float> q, vector<float> dq){
            vector<float> c;
            return c;
        }

        virtual float Ld_evaluate(vector<float> q, vector<float> fq){
            return 0.0;
        }


        virtual vector<float> solveD1Ld(vector<float> q, vector<float> p,vector<vector<particles>> sParticles){
            vector<float> fq;
            return fq;
        }

        virtual vector<float> solveD2Ld(vector<float> q,vector<float> fq,vector<vector<particles>> sParticles){
            vector<float> fp;
            return fp;
        }

};

class System1 : public L{
    public:
        int DOF = 3;
        vector<float> F = {2.0f, 2.0f, 0.0f};
        float k = 80.0f;
        float k1 = 80.0f;
        float m = 0.3f;
        float l = 5.0f;
        float h = 0.01f;

        System1(vector<float> q, vector<float> dq, vector<float> p,vector<vector<particles>> sParticles): L(q, dq, p) {

        }

        float L_evaluate(vector<float> q, vector<float> dq) override{
            return 0.5*m*(dq[0]*dq[0] + dq[1]*dq[1] + dq[2]*dq[2]) - 0.5*k1*pow((sqrt((q[0]*q[0] + q[1]*q[1] + q[2]*q[2])) - l) , 2) - (F[0]*q[0] + F[1]*q[1] + F[2]*q[2]);
        }     

        vector<float> dldq(vector<float> q, vector<float> dq) override{
            vector<float> c;
            for(int i=0;i<q.size();i++){
                c.push_back(0);
                c[i]= -k*q[i] - F[i];
            }
            return c;
        }

        vector<float> dldq_2(vector<float> q, vector<float> dq) override{
            vector<float> c;
            for(int i=0;i<q.size();i++){
                c.push_back(0);
                c[i]= m*dq[i];
            }
            return c;
        }

        float Ld_evaluate(vector<float> q, vector<float> fq) override{
            vector<float> t1;
            vector<float> t2;
            for(int i=0;i < q.size();i++){
                t1.push_back(0);
                t2.push_back(0);
            }
            for(int i=0;i < q.size();i++){
                t1[i]=(q[i]+fq[i])/2.0f;
                t2[i]=(fq[i]-q[i])/h;
            }
            return h*L_evaluate(t1,t2);
        }

        vector<float> solveD1Ld(vector<float> q, vector<float> p,vector<vector<particles>> sParticles) override{
            vector<float> fq;
            for(int k=0; k < q.size(); k++){
                fq.push_back(0);
            }
            fq[0] = ((((p[0]/h) - (F[0]*m/2))+((k1/2))+ q[0]*(4.0f*m - k*h*h) + 2.0f*k*h*h*l)/(4.0f*m+k*h*h));
            return fq;
        }

        vector<float> solveD2Ld(vector<float> q,vector<float> fq,vector<vector<particles>> sParticles) override{
            vector<float> fp;
            for(int i=0; i < q.size(); i++){
                fp.push_back(0);
            }
            fp[0] = -(((2.0f*h*h*F[0]*m)+ k*h*h*(q[0]+fq[0]-2.0f*l) - 4.0f*m*(-q[0]+fq[0]))/(4.0f*h));
            return fp;
        }        

};


class clothsim {
public:
    clothsim();

    void draw(void);
    void simulate(void);
    void applyeulermethod(void);
    void lagrangeintegrator(void);

    vector<vector<particles>> vParticles;
    vector<spring> vStrechsprings;
    vector<spring> vBendsprings;
};

clothsim::clothsim() {
    // particles
    for (int i = 0; i < 1; i++) {
        vector<particles> part;
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


void clothsim::draw(void) {
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


void clothsim::simulate(void){
    //Clear the force accumulator for all particles
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < MESH_SIZE; j++)
            vParticles[i][j].removeforces();
    }

    lagrangeintegrator();    
}


void clothsim::lagrangeintegrator(void){
    float h =TIME_STEP;
    float l = 5.0f;
    for (int i = 0; i < 1; i++) {
        for (int j = 1; j < MESH_SIZE; j++) {
            particles &p = vParticles[i][j];

            vector<float> q;
            q.push_back(p.m_pos.x);
            q.push_back(p.m_pos.y);
            q.push_back(p.m_pos.z);
            vector<float> dq;
            dq.push_back(p.m_vel.x);
            dq.push_back(p.m_vel.y);
            dq.push_back(p.m_vel.z);
            vector<float> pp;
            pp.push_back(p.m_mom.x);
            pp.push_back(p.m_mom.y);
            pp.push_back(p.m_mom.z);

            System1 sy(q,dq,pp,vParticles);

            L* sys = &sy;

            if(j==1){
                float a1 = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                float b1 = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                float c1 = sqrt(a1*a1 + b1*b1);
                float a2 = -(p.m_pos.x-vParticles[i][j+1].m_pos.x);
                float b2 = -(p.m_pos.y-vParticles[i][j+1].m_pos.y);
                float c2 = sqrt(a2*a2 + b2*b2);
                p.m_pos.x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].m_pos.x+vParticles[i][j+1].m_pos.x+l*(a1/c1)-l*(a2/c2)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].m_pos.x)+2.0f*(vParticles[i][j-1].m_pos.x)-(vParticles[i][j+2].m_pos.x)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].m_pos.y+vParticles[i][j+1].m_pos.y+l*(b1/c1)-l*(b2/c2)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].m_pos.y)+2.0f*(vParticles[i][j-1].m_pos.y)-(vParticles[i][j+2].m_pos.y)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.z = 0.0f;
                a1 = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                b1 = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                c1 = sqrt(a1*a1 + b1*b1);
                a2 = -(p.m_pos.x-vParticles[i][j+1].m_pos.x);
                b2 = -(p.m_pos.y-vParticles[i][j+1].m_pos.y);
                c2 = sqrt(a2*a2 + b2*b2);
                p.m_mom.x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.m_pos.x+q[0]-(vParticles[i][j-1].m_pos.x+vParticles[i][j+1].m_pos.x+l*(a1/c1)-l*(a2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.m_pos.x+q[0])+vParticles[i][j+2].m_pos.x-4.0f*(vParticles[i][j+1].m_pos.x)+2.0f*(vParticles[i][j-1].m_pos.x))-(mp/(h*h))*(p.m_pos.x-q[0]));
                p.m_mom.y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.m_pos.y+q[1]-(vParticles[i][j-1].m_pos.y+vParticles[i][j+1].m_pos.y+l*(b1/c1)-l*(b2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.m_pos.y+q[1])+vParticles[i][j+2].m_pos.y-4.0f*(vParticles[i][j+1].m_pos.y)+2.0f*(vParticles[i][j-1].m_pos.y))-(mp/(h*h))*(p.m_pos.y-q[1]));
                p.m_mom.z = 0.0f;
            }
            else if (j > 1 && j < MESH_SIZE - 2){
                float a1 = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                float b1 = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                float c1 = sqrt(a1*a1 + b1*b1);
                float a2 = -(p.m_pos.x-vParticles[i][j+1].m_pos.x);
                float b2 = -(p.m_pos.y-vParticles[i][j+1].m_pos.y);
                float c2 = sqrt(a2*a2 + b2*b2);
                p.m_pos.x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].m_pos.x+vParticles[i][j+1].m_pos.x+l*(a1/c1)-l*(a2/c2)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].m_pos.x+vParticles[i][j-1].m_pos.x)-(vParticles[i][j-2].m_pos.x+vParticles[i][j+2].m_pos.x)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].m_pos.y+vParticles[i][j+1].m_pos.y+l*(b1/c1)-l*(b2/c2)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(4.0f*(vParticles[i][j+1].m_pos.y+vParticles[i][j-1].m_pos.y)-(vParticles[i][j-2].m_pos.y+vParticles[i][j+2].m_pos.y)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(3.0f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.z = 0.0f;
                a1 = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                b1 = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                c1 = sqrt(a1*a1 + b1*b1);
                a2 = -(p.m_pos.x-vParticles[i][j+1].m_pos.x);
                b2 = -(p.m_pos.y-vParticles[i][j+1].m_pos.y);
                c2 = sqrt(a2*a2 + b2*b2);  
                p.m_mom.x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.m_pos.x+q[0]-(vParticles[i][j-1].m_pos.x+vParticles[i][j+1].m_pos.x+l*(a1/c1)-l*(a2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(3.0f*(p.m_pos.x+q[0])+vParticles[i][j-2].m_pos.x+vParticles[i][j+2].m_pos.x-4.0f*(vParticles[i][j+1].m_pos.x+vParticles[i][j-1].m_pos.x))-(mp/(h*h))*(p.m_pos.x-q[0]));
                p.m_mom.y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.m_pos.y+q[1]-(vParticles[i][j-1].m_pos.y+vParticles[i][j+1].m_pos.y+l*(b1/c1)-l*(b2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(3.0f*(p.m_pos.y+q[1])+vParticles[i][j-2].m_pos.y+vParticles[i][j+2].m_pos.y-4.0f*(vParticles[i][j+1].m_pos.y+vParticles[i][j-1].m_pos.y))-(mp/(h*h))*(p.m_pos.y-q[1]));
                p.m_mom.z = 0.0f;
            }
            else if(j==MESH_SIZE - 2){
                float a1 = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                float b1 = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                float c1 = sqrt(a1*a1 + b1*b1);
                float a2 = -(p.m_pos.x-vParticles[i][j+1].m_pos.x);
                float b2 = -(p.m_pos.y-vParticles[i][j+1].m_pos.y);
                float c2 = sqrt(a2*a2 + b2*b2);
                p.m_pos.x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].m_pos.x+vParticles[i][j+1].m_pos.x+l*(a1/c1)-l*(a2/c2)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j+1].m_pos.x)+4.0f*(vParticles[i][j-1].m_pos.x)-(vParticles[i][j-2].m_pos.x)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/2.0f)*(vParticles[i][j-1].m_pos.y+vParticles[i][j+1].m_pos.y+l*(b1/c1)-l*(b2/c2)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/2.0f)-(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j+1].m_pos.y)+4.0f*(vParticles[i][j-1].m_pos.y)-(vParticles[i][j-2].m_pos.y)))/((mp/(h*h)+(STRECH_SPRING_K/2.0f)+(2.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.z = 0.0f;
                a1 = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                b1 = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                c1 = sqrt(a1*a1 + b1*b1);
                a2 = -(p.m_pos.x-vParticles[i][j+1].m_pos.x);
                b2 = -(p.m_pos.y-vParticles[i][j+1].m_pos.y);
                c2 = sqrt(a2*a2 + b2*b2);
                p.m_mom.x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.m_pos.x+q[0]-(vParticles[i][j-1].m_pos.x+vParticles[i][j+1].m_pos.x+l*(a1/c1)-l*(a2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.m_pos.x+q[0])+vParticles[i][j-2].m_pos.x-2.0f*(vParticles[i][j+1].m_pos.x)-4.0f*(vParticles[i][j-1].m_pos.x))-(mp/(h*h))*(p.m_pos.x-q[0]));
                p.m_mom.y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/2.0f)*(p.m_pos.y+q[1]-(vParticles[i][j-1].m_pos.y+vParticles[i][j+1].m_pos.y+l*(b1/c1)-l*(b2/c2)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.5f*(p.m_pos.y+q[1])+vParticles[i][j-2].m_pos.y-2.0f*(vParticles[i][j+1].m_pos.y)-4.0f*(vParticles[i][j-1].m_pos.y))-(mp/(h*h))*(p.m_pos.y-q[1]));
                p.m_mom.z = 0.0f;
            }
            else if (j==MESH_SIZE - 1){
                float a = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                float b = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                float c = sqrt(a*a + b*b);  
                p.m_pos.x = ((((pp[0]/h)+(ACC_G[0]*mp/2.0f))+((STRECH_SPRING_K/4.0f)*(2.0f*vParticles[i][j-1].m_pos.x+2.0f*l*(a/c)))+ q[0]*((mp/(h*h))-(STRECH_SPRING_K/4.0f)-(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j-1].m_pos.x)-(vParticles[i][j-2].m_pos.x)))/((mp/(h*h)+(STRECH_SPRING_K/4.0f)+(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.y = ((((pp[1]/h)+(ACC_G[1]*mp/2.0f))+((STRECH_SPRING_K/4.0f)*(2.0f*vParticles[i][j-1].m_pos.y+2.0f*l*(b/c)))+ q[1]*((mp/(h*h))-(STRECH_SPRING_K/4.0f)-(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))+(BEND_SPRING_K/(2.0f*l*l*l))*(2.0f*(vParticles[i][j-1].m_pos.y)-(vParticles[i][j-2].m_pos.y)))/((mp/(h*h)+(STRECH_SPRING_K/4.0f)+(0.5f*BEND_SPRING_K/(2.0f*l*l*l)))));
                p.m_pos.z = 0.0f;

                a = (p.m_pos.x-vParticles[i][j-1].m_pos.x);
                b = (p.m_pos.y-vParticles[i][j-1].m_pos.y);
                c = sqrt(a*a + b*b);

                p.m_mom.x = -h*((-ACC_G[0]*mp/2.0f)+(STRECH_SPRING_K/4.0f)*(p.m_pos.x+q[0]-2.0f*(vParticles[i][j-1].m_pos.x+l*(a/c)))+(BEND_SPRING_K/(2.0f*l*l*l))*(0.5f*(p.m_pos.x+q[0])+vParticles[i][j-2].m_pos.x-2.0f*(vParticles[i][j-1].m_pos.x))-(mp/(h*h))*(p.m_pos.x-q[0]));
                p.m_mom.y = -h*((-ACC_G[1]*mp/2.0f)+(STRECH_SPRING_K/4.0f)*(p.m_pos.y+q[1]-2.0f*(vParticles[i][j-1].m_pos.y+l*(b/c)))+(BEND_SPRING_K/(2.0f*l*l*l))*(0.5f*(p.m_pos.y+q[1])+vParticles[i][j-2].m_pos.y-2.0f*(vParticles[i][j-1].m_pos.y))-(mp/(h*h))*(p.m_pos.y-q[1]));
                p.m_mom.z = 0.0f;

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

void SetTextures(const char *fileName) {
    
}
void SetLightSource(void) {
   
}

void SetMaterial(void) {
    
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(900, 700);
    glutCreateWindow("Cloth Simulation");

    // Initialize GLEW and other libraries

    SetTextures("cloth_texture.jpg");
    SetLightSource();
    SetMaterial();
    // Declare and initialize the cloth object

    glutDisplayFunc(DisplayCallback);
    glutIdleFunc([]() { cloth.simulate(); glutPostRedisplay(); });


    glutMainLoop();
    return 0;
}
