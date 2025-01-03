#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

# MIT License
# Copyright (c) 2022 Jakson-Almeida
# Consulte o arquivo LICENSE no repositório raiz para obter mais informações.

/* 16/09/2022 Biblioteca de Redes Neurais (C++, Arduino)
 * Autor: Almeida, Jakson Almeida
 * Contato: jakson.almeida@estudante.ufjf.br
 *
 * Redes Neurais Recorrentes
*/

class MatrizLin
{
  public:
    MatrizLin();
    MatrizLin(int mm, int nn);
    ~MatrizLin();
    void set(int i, int j, float val);
    float get(int i, int j);
    void printar();
    void printar(MatrizLin mat);
    int getNL();
    int getNC();
    void start();
  private:
    int nl, nc; // numero de linhas e colunas
    float *vet; // vetor de tamanho nl*nc
    int detInd(int i, int j);
};

MatrizLin::MatrizLin()
{
  nl = 1;
  nc = 1;
  vet = new float[nl*nc];
  start();
}

MatrizLin::MatrizLin(int mm, int nn)
{
  nl = mm;
  nc = nn;
  vet = new float[nl*nc];
  start();
}

MatrizLin::~MatrizLin()
{
  delete [] vet;
}

int MatrizLin::detInd(int i, int j)
{
  if(i >= 0 && i < nl && j >= 0 && j < nc)
    return i*nc + j;
  else
    return -1; // indice invalido
}

void MatrizLin::start()
{
  for(int i = 0; i < nl; i++)
    for(int j = 0; j < nc; j++)
      set(i, j, 0);
}

float MatrizLin::get(int i, int j)
{
  int k = detInd(i, j);
  if(k != -1)
    return vet[k];
  else {
    cout << "Indice invalido! (MatrizLin::get())" << endl;
  }
  return 0;
}

void MatrizLin::set(int i, int j, float val)
{
  int k = detInd(i, j);
  if(k != -1)
    vet[k] = val;
  else
    cout << "Indice invalido! (MatrizLin::set())" << endl;
}

void MatrizLin::printar()
{
  for(int i = 0; i < nl; i++) {
      cout << endl;
      for(int j = 0; j < nc; j++) {
        cout << get(i, j);
        cout << " \t";
      }
    }
    cout << endl;
}

void MatrizLin::printar(MatrizLin mat)
{
  int linhas = mat.getNL();
  int colunas = mat.getNC();
  for(int i = 0; i < linhas; i++) {
      cout << endl;
      for(int j = 0; j < colunas; j++) {
        cout << (mat.get(i, j));
        cout << " \t";
      }
    }
    cout << endl;
}

int MatrizLin::getNL()
{
  return nl;
}

int MatrizLin::getNC()
{
  return nc;
}

class Matrix {
  private:
    int linhas = 1;
    int colunas = 1;
  public:
    MatrizLin *mat;
    Matrix();
    Matrix(float m[2][2]);
    Matrix(float m[], int tm);
    Matrix(int n);
    Matrix(int l, int c);
    ~Matrix();
    void printMatrix();
    // void printMatrix2(float mat[2][2]);
    void printMatrix(Matrix A);
    Matrix matCopy();
    void transposta();
    void transposta(Matrix A);
    Matrix copy();
    Matrix copy(Matrix A);
    Matrix sum(Matrix A, Matrix B);
    void sum(Matrix A);
    Matrix sub(Matrix A, Matrix B);
    void sub(Matrix A);
    Matrix mult(Matrix A, float p);
    Matrix mult(Matrix A, Matrix B);
    Matrix* mult(Matrix *A, Matrix *B);
    void mult(Matrix *A, Matrix *B, Matrix *C);
    void mult(float p);
    void mult(Matrix A);
    void updateValues(float* values, int tm);
    void updateValues(int* values, int tm);
    void swapLines(int L0, int L1);
    void swapColumns(int c0, int c1);
    void moveDown(Matrix *mv);
    void moveDown2(Matrix mv);
    Matrix getLine(int li);
    Matrix getColunm(int col);
    int getNumLinhas();
    int getNumColunas();
    float get(int m, int n);
    bool pertence(int m, int n);
    void set(int m, int n, float v);
    void set(Matrix *mat);
    void setLine(Matrix lin, int li);
    void setColunm(Matrix colunm, int col);
};

Matrix::Matrix() {
  mat = new MatrizLin(linhas, colunas);
  mat->start();
}

Matrix::Matrix(float m[2][2]) {
  linhas = 2;
  colunas = 2;
  mat = new MatrizLin(linhas, colunas);
  for(int i = 0; i < linhas; i++)
    for(int j = 0; i < colunas; i++)
      mat->set(i, j, m[i][j]);
}

Matrix::Matrix(float m[], int tm) {
  linhas = 1;
  colunas = tm;
  mat = new MatrizLin(linhas, colunas);
  mat->start();
  cout << "Matriz vetor" << endl;
  for(int i = 0; i < colunas; i++) {
    mat->set(0, i, m[i]);
    cout << m[i] << endl;
  }
}

Matrix::Matrix(int n) {
  linhas = n;
  colunas = n;
  mat = new MatrizLin(linhas, colunas);
  mat->start();
}

Matrix::Matrix(int l, int c) {
  // cout << "Declaracao da Matriz" << endl;
  linhas = l;
  colunas = c;
  mat = new MatrizLin(linhas, colunas);
  mat->start();
}

Matrix::~Matrix() {
  mat->~MatrizLin();
}

void Matrix::printMatrix() {
  for(int i = 0; i < linhas; i++) {
      cout << endl;
      for(int j = 0; j < colunas; j++) {
        cout << get(i, j);
        cout << " \t";
      }
    }
    cout << endl;
}

// void Matrix::printMatrix2(float mat[2][2]) {
//   //Serial.println(linhas);
//   //Serial.println(colunas);
//   for(int i = 0; i < 2; i++) {
// //      Serial.println();
//       for(int j = 0; j < 2; j++) {
// //        Serial.print(mat[i][j]);
// //        Serial.print(" \t");
//       }
//     }
//     //Serial.println();
// }

void Matrix::printMatrix(Matrix A) {
  for(int i = 0; i < A.linhas; i++) {
      cout << endl;
      for(int j = 0; j < A.colunas; j++) {
        cout << A.mat->get(i, j);
        cout << " \t";
      }
    }
    cout << endl;
}

Matrix Matrix::matCopy() {
  Matrix C(linhas, colunas);
  for(int i = 0; i < linhas; i++) {
    for(int j = 0; j < colunas; j++) {
      C.mat->set(i, j, mat->get(i, j));
    }
  }
  return C;
}

void Matrix::transposta() {
  Matrix *C = new Matrix(colunas, linhas);
  // cout << "mat" << endl;
  // printMatrix();
  // cout << "antes" << endl;
  // C->printMatrix();
  for(int i = 0; i < linhas; i++) {
    for(int j = 0; j < colunas; j++) {
      C->mat->set(j, i, mat->get(i, j));
    }
  }
  mat = C->mat;
  // cout << "matriz" << endl;
  // C->printMatrix();
  int col = colunas;
  this->colunas = linhas;
  this->linhas = col;
}

Matrix Matrix::copy() {
  // cout << "copy" << endl;
  Matrix *C = new Matrix(linhas, colunas);
  // cout << "Matriz C" << endl;
  // cout << "Linhas: " << linhas << ", colunas: " << colunas << endl;
  for(int i = 0; i < linhas; i++) {
    for(int j = 0; j < colunas; j++) {
      C->mat->set(i, j, mat->get(i, j));
    }
  }
  // cout << "return" << endl;
  return *C;
}

Matrix Matrix::copy(Matrix A) {
  cout << "copy" << endl;
  Matrix *C = new Matrix(A.getNumLinhas(), A.getNumColunas());
  cout << "Matriz C" << endl;
  cout << "Linhas: " << A.getNumLinhas() << ", colunas: " << A.getNumColunas() << endl;
  for(int i = 0; i < A.getNumLinhas(); i++) {
    for(int j = 0; j < A.getNumColunas(); j++) {
      C->mat->set(i, j, A.mat->get(i, j));
    }
  }
  cout << "return" << endl;
  return *C;
}

Matrix Matrix::sum(Matrix A, Matrix B) {
  if((A.getNumColunas() != B.getNumColunas()) || (A.getNumLinhas() != B.getNumLinhas())) return A;
  Matrix *C = new Matrix(A.getNumColunas(), A.getNumLinhas());
  for(int i = 0; i < A.getNumLinhas(); i++) {
    for(int j = 0; j < A.getNumColunas(); j++) {
      C->mat->set(i, j, A.mat->get(i, j) + B.mat->get(i, j));
    }
  }
  return *C;
}

void Matrix::sum(Matrix A) {
  if((colunas != A.getNumColunas()) || (linhas != A.getNumLinhas())) {
    cout << "Erro ao somar" << endl;
    return;
  }
  // cout << "Somando" << endl;
  for(int i = 0; i < A.getNumLinhas(); i++) {
    // cout << endl;
    for(int j = 0; j < A.getNumColunas(); j++) {
      // cout << get(i, j);
      set(i, j, get(i, j) + A.get(i, j));
      // cout << " + ";
      // cout << A.get(i, j);
      // cout << " = ";
      // cout << get(i, j);
      // cout << endl;
    }
  }
  //cout << "Somou" << endl;
  // cout << endl;
}

Matrix Matrix::sub(Matrix A, Matrix B) {
  if((A.getNumColunas() != B.getNumColunas()) || (A.getNumLinhas() != B.getNumLinhas())) return A;
  Matrix *C = new Matrix(A.getNumColunas(), A.getNumLinhas());
  for(int i = 0; i < A.getNumLinhas(); i++) {
    for(int j = 0; j < A.getNumColunas(); j++) {
      C->mat->set(i, j, A.mat->get(i, j) - B.mat->get(i, j));
    }
  }
  return *C;
}

void Matrix::sub(Matrix A) {
  if((colunas != A.getNumColunas()) || (linhas != A.getNumLinhas())) return;
  for(int i = 0; i < A.getNumLinhas(); i++) {
    for(int j = 0; j < A.getNumColunas(); j++) {
      mat->set(i, j, mat->get(i, j) - A.mat->get(i, j));
    }
  }
}

Matrix Matrix::mult(Matrix A, float p) {
  Matrix *C = new Matrix(A.getNumColunas(), A.getNumLinhas());
  for(int i = 0; i < A.getNumLinhas(); i++) {
    for(int j = 0; j < A.getNumColunas(); j++) {
      C->mat->set(i, j, A.mat->get(i, j)*p);
    }
  }
  return *C;
}

Matrix Matrix::mult(Matrix A, Matrix B) {
  // cout << "Produto matricial" << endl;
  if((A.getNumColunas() != B.getNumLinhas())) {
    cout << "Erro em Matrix::mult(Matrix A, Matrix B)" << endl;
    return A;
  }
  float sum = 0;
  Matrix *C = new Matrix(A.getNumLinhas(), B.getNumColunas());
  for(int i = 0; i < A.getNumLinhas(); i++) {
    for(int j = 0; j < B.getNumColunas(); j++) {
      // float gambiarra = 0;
      sum = 0;
      // C->mat->set(i, j, 0);
      // cout << endl;
      for(int k = 0; k < A.getNumColunas(); k++) {
        // gambiarra = sum;
        sum += A.mat->get(i, k) * B.mat->get(k, j);
        // cout << gambiarra << " + " << A.mat->get(i, k) << " * " << B.mat->get(k, j) << " = " << sum << endl;
        // C->mat->set(i, j, C->mat->get(i, j) + A.mat->get(i, k) * B.mat->get(k, j));
      }
      C->mat->set(i, j, sum);
    }
  }
  // cout << endl << "O produto e igual a:" << endl;
  // C->printMatrix();
  return *C;
}

Matrix* Matrix::mult(Matrix *A, Matrix *B) {
  cout << "Produto matricial" << endl;
  if((A->getNumColunas() != B->getNumLinhas())) {
    cout << "Erro em Matrix::mult(Matrix A, Matrix B)" << endl;
    return A;
  }
  float sum = 0;
  Matrix *C = new Matrix(A->getNumLinhas(), B->getNumColunas());
  for(int i = 0; i < A->getNumLinhas(); i++) {
    for(int j = 0; j < B->getNumColunas(); j++) {
      float gambiarra = 0;
      sum = 0;
      // C->mat->set(i, j, 0);
      cout << endl;
      for(int k = 0; k < A->getNumColunas(); k++) {
        gambiarra = sum;
        sum += A->mat->get(i, k) * B->mat->get(k, j);
        cout << gambiarra << " + " << A->mat->get(i, k) << " * " << B->mat->get(k, j) << " = " << sum << endl;
        // C->mat->set(i, j, C->mat->get(i, j) + A.mat->get(i, k) * B.mat->get(k, j));
      }
      C->mat->set(i, j, sum);
    }
  }
  cout << endl << "O produto e igual a:" << endl;
  C->printMatrix();
  return C;
}

void Matrix::mult(Matrix *A, Matrix *B, Matrix *C) {
  // cout << "Produto matricial" << endl;
  if((A->getNumColunas() != B->getNumLinhas())) {
    cout << "Erro em Matrix::mult(Matrix A, Matrix B)" << endl;
    return;
  }
  C->~Matrix();
  C = new Matrix(A->getNumLinhas(), B->getNumColunas());
  float sum = 0;
  for(int i = 0; i < A->getNumLinhas(); i++) {
    for(int j = 0; j < B->getNumColunas(); j++) {
      // float gambiarra = 0;
      sum = 0;
      // C->mat->set(i, j, 0);
      // cout << endl;
      for(int k = 0; k < A->getNumColunas(); k++) {
        // gambiarra = sum;
        sum += A->mat->get(i, k) * B->mat->get(k, j);
        // cout << gambiarra << " + " << A.mat->get(i, k) << " * " << B.mat->get(k, j) << " = " << sum << endl;
        // C->mat->set(i, j, C->mat->get(i, j) + A.mat->get(i, k) * B.mat->get(k, j));
      }
      C->mat->set(i, j, sum);
    }
  }
  // cout << endl << "O produto e igual a:" << endl;
  // C->printMatrix();
  return;
}

void Matrix::mult(float p) {
  for(int i = 0; i < linhas; i++) {
    for(int j = 0; j < colunas; j++) {
      mat->set(i, j, mat->get(i, j)*p);
    }
  }
}

void Matrix::mult(Matrix A) {
  if((colunas != A.getNumLinhas())) return;
  Matrix *C = new Matrix(linhas, A.getNumColunas());
  for(int i = 0; i < linhas; i++) {
    for(int j = 0; j < A.getNumColunas(); j++) {
      C->mat->set(i, j, 0);
      for(int k = 0; k < colunas; k++) {
        C->mat->set(i, j, mat->get(i, k) + mat->get(i, k) * A.mat->get(k, j));
      }
    }
  }
  this->mat = C->mat;
  this->linhas = C->linhas;
  this->colunas = C->colunas;
}

void Matrix::updateValues(float* values, int tm) { //Update 06/04/2020
  if(this->linhas > 1 && this->colunas > 1) {
    cout << "Erro em updateValues" << endl;
    return;
  }
  int tamanho = tm;
  int minV = min(this->linhas, tamanho);
  if(this->linhas > this->colunas) {
    for(int i = 0; i < minV; i++) {
      this->mat->set(i, 0, values[i]);
    }
  }
  else {
    for(int i = 0; i < minV; i++) {
      this->mat->set(0, i, values[i]);
    }
  }
}

void Matrix::updateValues(int* values, int tm) { //Update 06/04/2020
  if(this->linhas > 1 && this->colunas > 1) return;
  int tamanho = tm;
  int minV = min(this->linhas, tamanho);
  if(this->linhas > this->colunas) {
    for(int i = 0; i < minV; i++) {
      this->mat->set(i, 0, values[i]);
    }
  }
  else {
    for(int i = 0; i < minV; i++) {
      this->mat->set(0, i, values[i]);
    }
  }
}

void Matrix::swapLines(int l0, int l1) {
  if(!pertence(l0, 0) || !pertence(l1, 0)) return;
  float li;
  for(int i = 0; i < colunas; i++) {
    li = mat->get(l0, i);
    mat->set(l0, i, mat->get(l1, i));
    mat->set(l1, i, li);
  }
}

Matrix Matrix::getLine(int li) {
  if(!pertence(li, 0)) {
    cout << "ERROR in Matrix::getLine(int li)" << endl;
    return Matrix();
  }
  Matrix *lin = new Matrix(1, colunas);
  for(int i = 0; i < lin->getNumColunas(); i++) lin->mat->set(0, i, this->mat->get(li, i));
  return *lin;
}

void Matrix::setLine(Matrix lin, int li) {
  if(getNumColunas() != lin.getNumColunas()) {
    cout << "ERROR in Matrix::setLine(Matrix lin, int li)";
    return;
  }
  if(li < 0 || li >= getNumLinhas()) {
    cout << "erro em if(li < 0 || li >= getNumLinhas())" << endl;
    return;
  }
  for(int i = 0; i < lin.getNumColunas(); i++) this->mat->set(li, i, lin.mat->get(0, i));
}

void Matrix::moveDown(Matrix *mv) {
  if(mv->getNumColunas() == colunas) {
    int ds = mv->getNumLinhas();
    if(ds > getNumLinhas()) {
      cout << "erro em Matrix::moveDown(Matrix mv)";
      return;
    }
    for(int i = linhas - 1; i >= ds; i--) {
      for(int j = 0; j < colunas; j++)
        set(i, j, get(i-ds, j));
    }
    for(int i = 0; i < ds; i++) {
      for(int j = 0; j < colunas; j++) {
        set(i, j, mv->get(i, j));
      }
    }
  }
  else cout << "erro em Matrix::moveDown(Matrix mv)";
}

void Matrix::moveDown2(Matrix mv) {
  if(mv.getNumColunas() == colunas) {
    int ds = mv.getNumLinhas();
    if(ds > getNumLinhas()) {
      cout << "erro em Matrix::moveDown(Matrix mv)";
      return;
    }
    for(int i = linhas - 1; i >= ds; i--) {
      for(int j = 0; j < colunas; j++)
        mat->set(i, j, mat->get(i-ds, j));
    }
    Matrix *m = new Matrix();
    for(int i = 0; i < ds; i++) {
      *m = mv.getLine(i);
      setLine(m->copy(), i);
    }
  }
  else cout << "erro em Matrix::moveDown(Matrix mv)";
}

Matrix Matrix::getColunm(int col) {
  if(!pertence(0, col)) return Matrix();
  Matrix *colunm = new Matrix(linhas, 1);
  for(int i = 0; i < colunm->getNumLinhas(); i++) colunm->mat->set(i, 0, mat->get(i, col));
  return *colunm;
}

int Matrix::getNumLinhas() {
  return this->linhas;
}

int Matrix::getNumColunas() {
  return this->colunas;
}

float Matrix::get(int m, int n) {
  return this->mat->get(m, n);
}

bool Matrix::pertence(int m, int n) {
  bool pert = (m < getNumLinhas() && n < getNumColunas() && m >= 0 && n >= 0);
  return pert;
}

void Matrix::set(int m, int n, float v) {
  this->mat->set(m, n, v);
}

void Matrix::set(Matrix *mat) {
  linhas = mat->getNumLinhas();
  colunas = mat->getNumColunas();
  this->mat = new MatrizLin(linhas, colunas);
  for(int i = 0; i < linhas; i++)
    for(int j = 0; j < colunas; j++)
      set(i, j, mat->get(i, j));
}

void Matrix::setColunm(Matrix colunm, int col) {
  if(getNumLinhas() != colunm.getNumLinhas()) return;
  if(col < 0 || col >= getNumColunas()) return;
  for(int i = 0; i < colunm.getNumLinhas(); i++) this->mat->set(i, col, colunm.mat->get(i, 0));
}

class Neural_Network {
  private:
    Matrix *layers;
    Matrix *deepLayers;
    Matrix *weights;
    Matrix vetValues;
    int *T_;
    int *presentLayers;
    int *tempLayers;
    int N_;
    int NP = 0;
    void start(int *values, int *tl, int TAM);
  public:
    int TAM = 0;
    Neural_Network();
    Neural_Network(int e, int s);
    Neural_Network(int e, int o, int s);
    Neural_Network(int values[], int tam);
    Neural_Network(int values[], int tl[], int tam);
    void debuga();
    void feedforward();
    void sigmoid(Matrix mat);
    void sigmoid(float vet[], int TAM);
    float sigmoid(float x);
    void reLu(Matrix *mat);
    void reLu(float vet[], int TAM);
    float reLu(float x);
    void updateValues(Matrix vet);
    void updateValues(float vet[], int tm);
    void updateEntrance(Matrix *entrance);
    void updateOutput(Matrix output);
    void updateHiddenWeights(Matrix weigths[]);
    Neural_Network get();
    void setVetValues(Matrix vet);
    int getEntranceLength();
    int getHiddenLength();
    int getOutputLength();
    Matrix getOutput();
    float* getOutputVet();
    float getOutput(int ind);
    Matrix getVetValues();
    int getVetValuesLength();
};

Neural_Network::Neural_Network() {
}

Neural_Network::Neural_Network(int e, int s) {
  int values[] = {e, s};
  int tl[2];
  start(values, tl, 2);
}

Neural_Network::Neural_Network(int e, int o, int s) {
  int values[] = {e, o, s};
  int tl[3];
  start(values, tl, 3);
}

Neural_Network::Neural_Network(int values[], int tam) {
  int *tl = new int[tam];
  for(int i = 0; i < tam; i++) tl[i] = 0;
  start(values, tl, tam);
}

Neural_Network::Neural_Network(int values[], int tl[], int tam) {
  start(values, tl, tam);
}

void Neural_Network::start(int values[], int tl[], int TAM) {
  this->TAM = TAM;
  int t1 = TAM;
  int t2 = TAM;
  // cout << "start" << endl;
//  Serial.println(TAM);
  if(t1 != t2 || t1 < 2) {
    cout << "ERRO no construtor Neural_Network(int[] values, int[] tl) e esperado vetores de tamanho iguais." << endl;
    return;
  }
  // cout << "aqui" << endl;
  this->presentLayers = new int[TAM];
  this->tempLayers    = new int[TAM];
  for(int i = 0; i < TAM; i++) {
    this->presentLayers[i] = values[i];
    this->tempLayers[i]    = tl[i];
  }
  // cout << "aqui" << endl;
  N_ = TAM;
  T_ = new int[N_];
  this->layers = new Matrix[N_];
  this->deepLayers = new Matrix[N_];
  this->weights = new Matrix[N_ - 1];
  // cout << "vetores" << endl;
  int v1, v2;
  Matrix *m4, *m3;
  for(int i = 0; i < N_; i++) {
    // cout << "for" << endl;
    T_[i] = presentLayers[i]*(tempLayers[i] + 1) + 1; // +1 => neur�nio adicional igual a 1 para o c�lculo dos bias
    // cout << 1 << endl;
    //Matrix *m3 = new Matrix((int)presentLayers[i], 1);
    // cout << 2 << endl;
    m4 = new Matrix(presentLayers[i], 1);
    m3 = new Matrix(T_[i], 1);
    layers[i] = *m4;
    // cout << 3 << endl;
    deepLayers[i] = *m3;
    // cout << 4 << endl;
  }
  // cout << "Aqui" << endl;
  for(int i = 0; i < N_ - 1; i++) {
    v1 = presentLayers[i+1];
    v2 = T_[i];
    Matrix *m2 = new Matrix(v1, v2);
    weights[i] = *m2;
    NP += v1*v2;
  }
  Matrix *m = new Matrix((int)NP, 1);
  this->vetValues = *m;
}

void Neural_Network::debuga() {
  cout << "DEBUGA REDE NEURAL" << endl;
  cout << endl;
  for(int i = 0; i < TAM; i++) {
    cout << presentLayers[i] << endl;
//    Serial.print(tempLayers[i]);
    cout << endl;
  }
  cout << "Print dos vetores LAYERS" << endl;
  for(int i = 0; i < N_; i++) {
    cout << i+1;
    cout << "_" << endl;
    // layers[i].printMatrix();
  }

//  Serial.println("Print dos vetores DEEPLAYERS");
//  for(int i = 0; i < N_; i++) {
//    Serial.print(i+1);
//    Serial.println("_");
//    deepLayers[i].printMatrix();
//  }

//  Serial.println("Print dos vetores WEIGHTS");
//  for(int i = 0; i < N_ - 1; i++) {
//    Serial.print(i+1);
//    Serial.println("_");
//    weights[i].printMatrix();
//  }
}

void Neural_Network::feedforward() {
  Matrix *mx = new Matrix();
  // layers[0].printMatrix(); cout << endl;
  for(int i = 0; i < N_ - 2; i++) {
    // cout << "Camada " + i << endl;
    //deepLayers[i].printMatrix();

    // cout << layers[i].get(0, 0);
    // cout << i << "_ numero de linhas " << layers[i].getNumLinhas() << " numero de colunas " << layers[i].getNumColunas() << endl;

    cout << "moveDown(&layers[" << i << "])" << endl;
    deepLayers[i].moveDown(&layers[i]);
    cout << "setou final igual a 1" << endl;
    deepLayers[i].set(deepLayers[i].getNumLinhas() - 1, 0, 1);

    cout << "Produto das duas matrizes" << endl;
    weights[i].printMatrix(); cout << endl;
    deepLayers[i].printMatrix(); cout << endl;

    // layers[i+1].~Matrix();
    mx = layers[i+1].mult(&weights[i], &deepLayers[i]);
    cout << "print mx" << endl;
    mx->printMatrix();
    layers[i+1] = *mx;
    // layers[i+1].mult(&weights[i], &deepLayers[i], &layers[i+1]);

    cout << "o produto e:" << endl;
    layers[i+1].printMatrix(); cout << endl;
    // cout << "Mat1 " << weights[i].getNumLinhas() << " " << weights[i].getNumColunas() <<endl;
    // cout << "Mat2 " << deepLayers[i].getNumLinhas() << " " << deepLayers[i].getNumColunas() <<endl;
    reLu(&layers[i+1]);

    // cout << "reLu" << endl;
    // layers[i+1].printMatrix(); cout << endl;
    // sigmoid(layers[i+1]);
    //deepLayers[i].printMatrix(); println();
  }
  int ind  = N_ - 2;
  int ind2 = N_ - 1;

  cout << "moveDown(&layers[" << ind << "])" << endl;
  deepLayers[ind].moveDown(&layers[ind]);
  cout << "setou final igual a 1" << endl;
  deepLayers[ind].set(deepLayers[ind].getNumLinhas() - 1, 0, 1);

  cout << "Produto das duas matrizes" << endl;
  weights[ind].printMatrix(); cout << endl;
  deepLayers[ind].printMatrix(); cout << endl;

  // layers[ind2].~Matrix();
  mx = layers[ind2].mult(&weights[ind], &deepLayers[ind]);
  cout << "print mx" << endl;
  mx->printMatrix();
  layers[ind2] = *mx;
  // layers[ind2].mult(&weights[ind], &deepLayers[ind], &layers[ind2]);
  cout << "o produto e:" << endl;
  layers[ind2].printMatrix();
}

void Neural_Network::sigmoid(Matrix mat) {
  float *vet;
  int col = mat.getNumColunas();
  int lin = mat.getNumLinhas();
  if (col > 0) {
    vet = new float[col];
    for(int i = 0; i < col; i++) {
      vet[i] = mat.mat->get(0, i);
      vet[i] = sigmoid(vet[i]);
    }
  }
else {
    vet = new float[lin];
    mat.transposta();
    for(int i = 0; i < col; i++) {
      vet[i] = mat.mat->get(0, i);
      vet[i] = sigmoid(vet[i]);
    }
  }
  mat.transposta();
}

void Neural_Network::sigmoid(float vet[], int TAM) {
  for (int i = 0; i < TAM; i++) {
    vet[i] = sigmoid(vet[i]);
  }
}

float Neural_Network::sigmoid(float x) {
  return 1.0 / (1.0 + exp(-x));
}

void Neural_Network::reLu(Matrix *mat) {
  for(int i = 0; i < mat->getNumLinhas(); i++) {
    for(int j = 0; j < mat->getNumColunas(); j++)
      mat->set(i, j, reLu(mat->get(i, j)));
  }
}

void Neural_Network::reLu(float vet[], int TAM) {
  for (int i = 0; i < TAM; i++) {
    vet[i] = reLu(vet[i]);
  }
}

float Neural_Network::reLu(float x) {
  return max((float) 0.0, x);
}

void Neural_Network::updateValues(Matrix vet) {
  vetValues = vet.copy();
  if(vet.getNumLinhas() != NP) {
    cout << "Error in updateValues" << endl;
    return;
  }
  int ind = 0;
  for(int i = 0; i < N_-1; i++) {
    for(int m = 0; m < weights[i].getNumLinhas(); m++) {
      for(int n = 0; n < weights[i].getNumColunas(); n++) {
        weights[i].set(m, n, vet.get(ind, 0));
        ind++;
      }
    }
  }
}

void Neural_Network::updateValues(float vet[], int tm) {
  Matrix *mat = new Matrix(vet, tm);
  mat->transposta();
  updateValues(*mat);
}

void Neural_Network::updateEntrance(Matrix *entrance) {
  this->layers[0] = *entrance;
}

void Neural_Network::updateOutput(Matrix output) {
  this->layers[N_-1] = output;
}

void Neural_Network::updateHiddenWeights(Matrix weights[]) {
  this->weights = weights;
}

Neural_Network Neural_Network::get() {
  Neural_Network net(presentLayers, tempLayers, TAM);
  net.updateValues(vetValues.copy());
  return net;
}

void Neural_Network::setVetValues(Matrix vet) {
  if (vet.getNumLinhas() == vetValues.getNumLinhas() && vet.getNumColunas() == 1) this->vetValues = vet;
  else cout << "Error in setVetValues()" << endl;
}

int Neural_Network::getEntranceLength() {
  return this->layers[0].getNumLinhas();
}

int Neural_Network::getHiddenLength() {
  return N_-2;
}

int Neural_Network::getOutputLength() {
  return this->layers[N_-1].getNumLinhas();
}

Matrix Neural_Network::getOutput() {
  return this->layers[N_-1];
}

float* Neural_Network::getOutputVet() {
  float *vet = new float[layers[N_-1].getNumLinhas()];
  for (int i = 0; i < layers[N_-1].getNumLinhas(); i++) {
    vet[i] = layers[N_-1].mat->get(i, 0);
  }
  return vet;
}

float Neural_Network::getOutput(int ind) {
  return layers[N_-1].mat->get(ind, 0);
}

int Neural_Network::getVetValuesLength() {
  return this->vetValues.getNumLinhas();
}

void printSum(Matrix A) {
  Matrix *C = new Matrix(A.getNumLinhas(), A.getNumColunas());
  cout << "Somando e debugando" << endl;
  for(int i = 0; i < A.getNumLinhas(); i++) {
    cout << endl;
    for(int j = 0; j < A.getNumColunas(); j++) {
      cout << C->get(i, j);
      C->set(i, j, C->get(i, j) + A.get(i, j));
      cout << " + ";
      cout << A.get(i, j);
      cout << " = ";
      cout << C->get(i, j);
      cout << endl;
    }
  }
  //cout << "Somou" << endl;
  cout << endl;
  cout << "Bora printar?" << endl;
  cout << "So se for agora:" << endl;
  C->printMatrix();
  cout << "(Tudo igual)" << endl;
  A.printMatrix();
  cout << "(Sem mudanca)" << endl;
}

float *out;
int deepLearning[]   = {2, 10, 10, 2};
int temporalLayers[] = {2, 1, 0, 0};
float vet[] = {1, -1};

int main()
{
    int TAM = sizeof(deepLearning) / sizeof(int);
    int tamVet = sizeof(vet) / sizeof(float);
    Neural_Network neuron(deepLearning, temporalLayers, TAM);
    Matrix weightsAndBias(neuron.getVetValuesLength(), 1);
    Matrix vetEntrance(tamVet, 1);
    for(int i = 0; i < weightsAndBias.getNumLinhas(); i++)
      weightsAndBias.set(i, 0, i*0.3 - 3);
    neuron.updateValues(weightsAndBias);
    vetEntrance.updateValues(vet, tamVet);
    neuron.updateEntrance(&vetEntrance);
    neuron.feedforward();
    out = neuron.getOutputVet();
    // cout << "Primeiro" << endl;

    std::cout << std::fixed;
    std::cout << std::setprecision(2);

    cout << "v1 = " << out[0] << endl;
    cout << "v2 = " << out[1] << endl << endl;

    cout << "Segundo" << endl;
    vet[0] = 0;
    vet[1] = 0;
    vetEntrance.updateValues(vet, tamVet);
    neuron.updateEntrance(&vetEntrance);
    neuron.feedforward();
    cout << "v1 = " << neuron.getOutput(0) << endl;
    cout << "v2 = " << neuron.getOutput(1) << endl;

    return 0;
}
