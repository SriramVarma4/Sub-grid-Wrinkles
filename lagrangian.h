#include <vector>
#include <GL/glut.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>


#ifndef LAGRANGIAN_H

#define LAGRANGIAN_H

class particles {
    private:
        const float mass;
        glm::vec3 position;
        glm::vec3 velocity{0.0f, 0.0f, 0.0f};
        glm::vec3 momentum{0.0f, 0.0f, 0.0f};
        glm::vec3 acceleration{0.0f, 0.0f, 0.0f};
        glm::vec3 force{0.0f, 0.0f, 0.0f};
        glm::vec2 texture;
        glm::vec3 normal{0.0f, 0.0f, 0.0f};

    public:
        particles(float m, glm::vec3 pos, glm::vec2 tex) : mass{m}, position{pos}, texture{tex}{

        }

        glm::vec3 const& getposition() const{
            return position;
        }

        glm::vec3 const& getvelocity() const{
            return velocity;
        }

        glm::vec3 const& getmomentum() const{
            return momentum;
        }

        void draw(){
            glColor3f(1.0f, 1.0f, 1.0f);
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glVertex3fv(glm::value_ptr(position));
            glEnd();
        }

        void removeforces(){
            force = glm::vec3(0.0f, 0.0f, 0.0f);
        }

        void setposition(const glm::vec3& p) {
           position = p;
        }

        void setvelocity(const glm::vec3& v) {
            velocity = v;
        }

        void setmomentum(const glm::vec3& m) {
            momentum = m;
        }

};

class Lagrangian{
    protected:
        glm::vec3 q_k;
        glm::vec3 dq_k;
        glm::vec3 p_k;

    public:
        Lagrangian(glm::vec3 q, glm::vec3 dq, glm::vec3 p) : q_k{q}, dq_k{dq}, p_k{p}{

        }

        virtual ~Lagrangian(){

        }

        virtual float L_evaluate(const glm::vec3& q, const glm::vec3& dq) = 0;

        virtual glm::vec3 dldq(const glm::vec3& q, const glm::vec3& dq) = 0;

        virtual glm::vec3 dldq_2(const glm::vec3& q, const glm::vec3& dq) = 0;

        // virtual float Ld_evaluate(const glm::vec3& q, const glm::vec3& fq) = 0;

        virtual std::vector<float> solveD1Ld(const std::vector<float>& q, const std::vector<float>& p, const std::vector<std::vector<particles>>& sParticles) = 0;

        virtual std::vector<float> solveD2Ld(const std::vector<float>& q, const std::vector<float>& fq, const std::vector<std::vector<particles>>& sParticles) = 0;

        virtual void run() = 0;
};


#endif