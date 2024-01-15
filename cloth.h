#include "lagrangian.h"

#ifndef CLOTH_H

#define CLOTH_H

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
        glVertex3fv(glm::value_ptr(m_particle1->getposition()));
        glVertex3fv(glm::value_ptr(m_particle2->getposition()));
        glEnd();
    }

    float deflength(){
        return glm::length(m_particle1->getposition() - m_particle2->getposition());
    }

    glm::vec3 springforce(){
        glm::vec3 link = m_particle2->getposition() - m_particle1->getposition();
        float slength = glm::length(link);
        glm::vec3 force = link*(1 - m_L/slength) * m_K;
        return force;
    }
};

class clothsim {
public:
    clothsim();

    void draw();

    void simulate();

    void applyeulermethod();

    void lagrangeintegrator();

    std::vector<std::vector<particles>> vParticles;
    std::vector<spring> vStrechsprings;
    std::vector<spring> vBendsprings;
};

#endif