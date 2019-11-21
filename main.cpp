#include <iostream>
#include <cmath>
#include <stdio.h>
#include <string>

using namespace std;

const int lengthX = 11;
const int lengthY = 60;
const int gap = 9;

const int q = 9;
const float w1 = 4.0 / 9.0;
const float w2 = 1.0 / 9.0;
const float w3 = 1.0 / 36.0;
const float c_squ = 1.0 / 3.0;
const float deltaU = 1.32e-3;
const float omega = 1.0;
float u = 0.08;
const int index_map[] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

float rho[lengthX][lengthY];
float f[lengthX][lengthY][q];
int bound[lengthX][lengthY];
int numActiveNodes = 0;

float computeDisFunEq(float, float, float, float);

void setupColorsMap();

void shiftRightX(int);

void shiftRightY(int);

void shiftLiftX(int);

void shiftLiftY(int);

void lbm(bool, float);

void puazeyl(float, float, int);

void tube();

void barrier();

float funVelocity(float, float);

float uMax = 0.0f;
uint32_t color[1024];
//string color[1024];

/*
 *
 */
int main(int argc, char **argv) {
    setupColorsMap();
    float st;

    //numActiveNodes=0;
    //tube();
    //lbm(false,1);
    //st = (uMax/1024);
    //printf("\n Max :%.3f \n st:%f\n\n",uMax,st);
    //numActiveNodes=0;
    //tube();
    //lbm(true,st);
    //printf("\n Max :%.3f \n st:%f\n\n",uMax,st);createColor
    //puazeyl(deltaU, omega, axisX);

    numActiveNodes = 0;
    barrier();
    uMax = 0.0f;
    lbm(false, 0.1f);
    st = uMax / 1024.0f;
    numActiveNodes = 0;
    barrier();
    uMax = 0.0f;
    lbm(true, st);

    return 0;
}

void lbm(bool statusPrint, float st) {
    float velocityX[lengthX][lengthY];
    float velocityY[lengthX][lengthY];
    float u_prev = 0;
    int t = 0;

    while ((t < 10000000 && abs(u - u_prev) >= 1e-7) || t < 100) {

        shiftRightX(1);
        shiftRightY(2);
        shiftLiftX(3);
        shiftLiftY(4);
        shiftRightX(5);
        shiftRightY(5);
        shiftLiftX(6);
        shiftRightY(6);
        shiftLiftX(7);
        shiftLiftY(7);
        shiftRightX(8);
        shiftLiftY(8);

        float tmp[lengthX][lengthY][q];
        for (int i = 0; i < lengthX; i++)
            for (int j = 0; j < lengthY; j++) {
                rho[i][j] = 0;
                for (int k = 0; k < q; k++) {
                    tmp[i][j][k] = f[i][j][k];
                    rho[i][j] += f[i][j][k];
                }
            }

        for (int i = 0; i < lengthX; i++)
            for (int j = 0; j < lengthY; j++) {
                float xR = 0, xL = 0, yR = 0, yL = 0;
                for (int k : {1, 5, 8}) xR += f[i][j][k];
                for (int k : {3, 6, 7}) xL += f[i][j][k];
                for (int k : {2, 5, 6}) yR += f[i][j][k];
                for (int k : {4, 7, 8}) yL += f[i][j][k];
                velocityX[i][j] = (xR - xL) / rho[i][j];
                velocityY[i][j] = (yR - yL) / rho[i][j] + deltaU;
            }

        for (int i = 0; i < lengthX; i++)
            for (int j = 0; j < lengthY; j++) {
                if (bound[i][j] == 1) velocityX[i][j] = 0.0;
                if (bound[i][j] == 1) velocityY[i][j] = 0.0;
                if (bound[i][j] == 1) rho[i][j] = 0.0;

                float u_squ, u_c2, u_c4, u_c6, u_c8;
                u_squ = velocityX[i][j] * velocityX[i][j] + velocityY[i][j] * velocityY[i][j];
                u_c2 = velocityX[i][j] + velocityY[i][j];
                u_c4 = -velocityX[i][j] + velocityY[i][j];
                u_c6 = -u_c2;
                u_c8 = -u_c4;
                float feq[q];
                feq[0] = computeDisFunEq(w1, rho[i][j], 0.0, u_squ);
                feq[1] = computeDisFunEq(w2, rho[i][j], velocityX[i][j], u_squ);
                feq[2] = computeDisFunEq(w2, rho[i][j], velocityY[i][j], u_squ);
                feq[3] = computeDisFunEq(w2, rho[i][j], -velocityX[i][j], u_squ);
                feq[4] = computeDisFunEq(w2, rho[i][j], -velocityY[i][j], u_squ);
                feq[5] = computeDisFunEq(w3, rho[i][j], u_c2, u_squ);
                feq[6] = computeDisFunEq(w3, rho[i][j], u_c4, u_squ);
                feq[7] = computeDisFunEq(w3, rho[i][j], u_c6, u_squ);
                feq[8] = computeDisFunEq(w3, rho[i][j], u_c8, u_squ);
                for (int k = 0; k < q; k++)
                    f[i][j][k] = omega * feq[k] + (1 - omega) * f[i][j][k];
            }

        for (int i = 0; i < lengthX; i++)
            for (int j = 0; j < lengthY; j++) {
                if (bound[i][j] == 0) continue;
                for (int k = 0; k < q; k++) {
                    f[i][j][k] = tmp[i][j][index_map[k]];
                }
            }
        u_prev = u;
        u = 0;
        for (int i = 0; i < lengthX; i++)
            for (int j = 0; j < lengthY; j++)
                u += sqrt(velocityY[i][j] * velocityY[i][j] + velocityX[i][j] * velocityX[i][j]);
        u = u / numActiveNodes;
        t++;
        for (int i = 0; i < lengthY; i++) {
            for (int j = 0; j < lengthX; j++) {
                float v = funVelocity(velocityX[j][i], velocityY[j][i]);
                if (v > uMax) uMax = v;
                if (statusPrint) {
                    int c = abs((int) (v / st) - 1);
                    if (c > 1023) c = 1023;
                    printf(" %d %d %d", j, i, (bound[j][i] == 1) ? 0 : color[c]);
                    //printf(" %d %d %s",j,i,(bound[j][i] == 1) ? "0 0 0": color[c].c_str());
                    //printf(" %f % 4d", v,c);
                }

            }
        }
        if (statusPrint) printf("\n");
    }
    //printf("\nMax : %f t:%d st :%f\n\n",uMax,t,st);
}

void puazeyl(float delta, float omegag, int width) {

    float g = delta * omega;
    float nu = 1.0f / 6.0f;
    float a = ((float) width - 2.0f) / 2.f;
    int k = width / 2;
    float arr[width];
    for (int i = 0; i <= k; i++) {
        if (k + i >= width - 1) {
            arr[k + i] = a;
            arr[k - i] = -a;
        } else {
            arr[k + i] = (float) i;
            arr[k - i] = (float) -i;
        }
    }
    printf("\n");
    for (float x : arr) {
        printf(" %.3f", g * (a * a - x * x) / (2.0f * nu));
    }
    printf("\n");
}

float funVelocity(float x, float y) {
    float us = sqrt(x * x + y * y);
    return us;
}

float computeDisFunEq(float w, float rho, float u, float u_squ) {
    float m = u / c_squ;
    float l = w * rho * (1.0f + m + 0.5f * m * m - u_squ / (2.0f * c_squ));
    return l;
}

void shiftRightX(const int k) {
    for (int j = 0; j < lengthY; j++) {
        float t = f[lengthX - 1][j][k];
        for (int i = lengthX - 1; i >= 1; i--)
            f[i][j][k] = f[i - 1][j][k];
        f[0][j][k] = t;
    }
}

void shiftRightY(const int k) {
    for (int i = 0; i < lengthX; i++) {
        float t = f[i][lengthY - 1][k];
        for (int j = lengthY - 1; j >= 1; j--)
            f[i][j][k] = f[i][j - 1][k];
        f[i][0][k] = t;
    }
}

void shiftLiftX(const int k) {
    for (int j = 0; j < lengthY; j++) {
        float t = f[lengthX - 1][j][k];
        for (int i = 0; i < lengthX - 1; i++)
            f[i][j][k] = f[i + 1][j][k];
        f[lengthX - 1][j][k] = t;
    }
}

void shiftLiftY(const int k) {
    for (int i = 0; i < lengthX; i++) {
        float t = f[i][lengthX - 1][k];
        for (int j = 0; j < lengthY - 1; j++)
            f[i][j][k] = f[i][j + 1][k];
        f[i][lengthY - 1][k] = t;
    }
}

void tube() {
    for (int i = 0; i < lengthX; i++) {
        for (int j = 0; j < lengthY; j++) {
            bound[i][j] = 0;
            if (i == 0 || i == lengthX - 1) bound[i][j] = 1;
            rho[i][j] = 1.0f;
            numActiveNodes += 1 - bound[i][j];
            for (int k = 0; k < q; k++) {
                f[i][j][k] = rho[i][j] / 9.0f;
            }
        }
    }
}

void barrier() {
    int r = lengthX - 1;

    for (int i = 0; i < lengthX; i++)
        for (int j = 0; j < lengthY - lengthX - gap; j++) {
            bound[i][j] = 0;
            if (i == 0 || i == lengthX - 1) bound[i][j] = 1;
            numActiveNodes += 1 - bound[i][j];
        }

    for (int i = 0; i < lengthX; i++)
        for (int j = lengthY - lengthX - gap; j < lengthY - lengthX; j++) {
            bound[i][j] = 0;
            numActiveNodes += 1 - bound[i][j];
        }

    for (int i = 0; i < lengthX; i++)
        for (int j = lengthY - lengthX; j < lengthY; j++) {
            bound[i][j] = 0;
            if ((i * i + (lengthY - j - 1) * (lengthY - j - 1)) <= r * r) bound[i][j] = 1;
            numActiveNodes += 1 - bound[i][j];
        }

    for (int i = 0; i < lengthX; i++)
        for (int j = 0; j < lengthY; j++) {
            for (int k = 0; k < q; k++) {
                rho[i][j] = 1.0f;
                f[i][j][k] = rho[i][j] / 9.0f;
            }
        }
}

void setupColorsMap() {
    int k = 0;
    for (int i = 0; i <= 255; i++) {
        //char s[11];sprintf(s, "%d %d %d", 0, 255, i);color[k++] = string(s);
        color[k++] = 0 << 8 * 2 | 255 << 8 | i;
    }
    for (int i = 255; i >= 0; i--) {
        //char s[11];sprintf(s, "%d %d %d", 0, i, 255);color[k++] = string(s);
        color[k++] = 0 << 8 * 2 | i << 8 | 255;
    }
    for (int i = 0; i <= 255; i++) {
        //char s[11];sprintf(s, "%d %d %d", i, 0, 255);color[k++] = string(s);
        color[k++] = i << 8 * 2 | 0 << 8 | 255;
    }
    for (int i = 255; i >= 0; i--) {
        //char s[11];sprintf(s, "%d %d %d", 255, 0, i);color[k++] = string(s);
        color[k++] = 255 << 8 * 2 | 0 << 8 | i;
    }
}
