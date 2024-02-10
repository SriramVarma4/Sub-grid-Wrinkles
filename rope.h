#include <GL/glut.h>
#include <glm/glm.hpp>
#include <vector>
#include <glm/gtc/type_ptr.hpp>

#ifndef CLOTH_H

#define CLOTH_H


#define mp 0.3f
#define STRECH_SPRING_L 5.0f
#define STRECH_SPRING_K 80.0f
#define ROPE_SIZE 20

#define ACC_G glm::vec3(-2.0f,-2.0f,0.0f)
#define TIME_STEP 0.01f

class particles {
    private:
        const float mass;

        // Particle properties
        glm::vec3 pos;
        glm::vec3 vel{0.0f, 0.0f, 0.0f};
        glm::vec3 mom{0.0f, 0.0f, 0.0f};

        // Derived properties
        glm::vec3 acceleration{0.0f, 0.0f, 0.0f};
        glm::vec3 force{0.0f, 0.0f, 0.0f};

        glm::vec2 texture;
        glm::vec3 normal{0.0f, 0.0f, 0.0f};

    public:
        particles(float m, glm::vec3 pos, glm::vec2 tex) : mass{m}, pos{pos}, texture{tex}{

        }

        glm::vec3& position(){
            return pos;
        }

        glm::vec3& velocity(){
            return vel;
        }

        glm::vec3& momentum(){
            return mom;
        }

        void draw(){
            glColor3f(1.0f, 1.0f, 1.0f);
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glVertex3fv(glm::value_ptr(pos));
            glEnd();
        }

        void removeforces(){
            force = glm::vec3(0.0f, 0.0f, 0.0f);
        }

};

class spring {
private:
    particles *m_particle1;
    particles *m_particle2;
    const float m_L;
    const float m_K;
public:
    spring(particles *p1, particles *p2, float L, float K) : m_L{L}, m_K{K}, m_particle1{p1}, m_particle2{p2}{

    }

    void draw(){
        glColor3f(1.0f, 1.0f, 1.0f);
        glBegin(GL_LINES);
        glVertex3fv(glm::value_ptr(m_particle1->position()));
        glVertex3fv(glm::value_ptr(m_particle2->position()));
        glEnd();
    }

    float deflength(){
        return glm::length(m_particle1->position() - m_particle2->position());
    }

    glm::vec3 springforce(){
        glm::vec3 link = m_particle2->position() - m_particle1->position();
        float slength = glm::length(link);
        glm::vec3 force = link*(1 - m_L/slength) * m_K;
        return force;
    }
};

class rope{
public:
    // We have particles in cloth
    std::vector<particles> vParticles;

    // We have stretch springs in cloth
    std::vector<spring> vStrechsprings;

    // We have bend springs in cloth
    std::vector<spring> vBendsprings;

    rope(){
        glm::vec3 origin(0.0f, 0.0f, 0.0f);

        // Initial position of all particles on rope is on x axis seperated by their natural length of springs
        for (int j = 0; j < ROPE_SIZE; j++) {

            glm::vec2 texcoord(j / (ROPE_SIZE - 1.0f), 1 / (ROPE_SIZE - 1.0f));
            //glm::vec2 texcoord(0.0f, 0.0f);
            vParticles.push_back(particles(mp, origin, texcoord));
            origin += glm::vec3(STRECH_SPRING_L, 0.0f, 0.0f);

        }

        for (int j = 0; j < ROPE_SIZE - 1; j++){
            vStrechsprings.push_back(spring(&vParticles[j], &vParticles[j + 1], STRECH_SPRING_L, STRECH_SPRING_K));
        }
    }

    void draw(){
        // Clear the color buffer and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw springs
        for (int i = 0; i < vStrechsprings.size(); i++) {
            vStrechsprings[i].draw();
        }

        // Draw particles
        for (int j = 0; j < ROPE_SIZE; j++) {
            vParticles[j].draw();
        }
    
    }

    void simulate(){
        //Clear the force accumulator for all particles
        for (int j = 0; j < ROPE_SIZE; j++){
            vParticles[j].removeforces();
        }

        timeintegrator();
    }


    void timeintegrator(){


        for (int j = 1; j < ROPE_SIZE; j++) {

            particles& p = vParticles[j];

            glm::vec3& q = p.position();
            glm::vec3& q_dot = p.velocity();
            glm::vec3 momentum = p.momentum();

        }


    }

};

#endif