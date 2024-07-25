#include <iostream>
#include <corecrt_math_defines.h>
#include <vector>

using namespace std;
using Matrix3x3 = vector<vector<double>>;
// 例えば以下のようにMatrix3x3型の変数を定義することができる
// Matrix3x3 mat = {
//     {1, 2, 3},
//     {0, 1, 4},
//     {5, 6, 0}
// };
// 科学計算なので全部double型で計算😎

// 行列の掛け算を行う関数
void multiplyMatrix(const Matrix3x3& matrixA, const Matrix3x3& matrixB, Matrix3x3& result) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = 0.0f;
            for (int k = 0; k < 3; ++k) {
                result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
}

// 出力のための関数
void printMatrix(const Matrix3x3& matrix) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
void printVector(const vector<double>& vector) {
    for (int i = 0; i < 3; ++i) {
        cout << vector[i] << " ";
    }
        cout << endl;
}

// 2つのベクトルの内積を計算する関数
double dotProduct(const vector<double>& vec1, const vector<double>& vec2) {
    double dot = 0.0;
    for (size_t i = 0; i < vec1.size(); ++i) {
        dot += vec1[i] * vec2[i];
    }
    return dot;
}

// ベクトルの大きさを計算する関数
double magnitude(const vector<double>& vec) {
    double mag = 0.0;
    for (double val : vec) {
        mag += val * val;
    }
    return sqrt(mag);
}

// 2つのベクトルのなす角を計算する関数
double angleBetweenVectors(const vector<double>& vec1, const vector<double>& vec2) {
    double dot = dotProduct(vec1, vec2);
    double mag1 = magnitude(vec1);
    double mag2 = magnitude(vec2);
    double cosTheta = dot / (mag1 * mag2);
    return acos(cosTheta) * (180.0 / M_PI);
}

// 3x3行列の行列式を計算する関数
double determinant(const Matrix3x3& mat) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
        - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
        + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

// 3x3行列の余因子行列を計算する関数
Matrix3x3 adjoint(const Matrix3x3& mat) {
    Matrix3x3 adj(3, std::vector<double>(3));

    adj[0][0] = mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1];
    adj[0][1] = mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2];
    adj[0][2] = mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1];

    adj[1][0] = mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2];
    adj[1][1] = mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0];
    adj[1][2] = mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2];

    adj[2][0] = mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0];
    adj[2][1] = mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1];
    adj[2][2] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

    return adj;
}

// 3x3行列の逆行列を計算する関数
Matrix3x3 inverse(const Matrix3x3& mat) {
    double det = determinant(mat);
    if (det == 0) {
        throw runtime_error("Matrix is singular and cannot be inverted.");
    }

    Matrix3x3 adj = adjoint(mat);
    Matrix3x3 inv(3, vector<double>(3));

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            inv[i][j] = adj[i][j] / det;
        }
    }

    return inv;
}



// 24個の固定行列を定義
const vector<Matrix3x3> fixedMatrices = {
    {{1, 0, 0},  {0, 1, 0},  {0, 0, 1}},  // V1
    {{1, 0, 0},  {0, -1, 0}, {0, 0, -1}}, // V2
    {{-1, 0, 0}, {0, 1, 0},  {0, 0, -1}}, // V3
    {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},  // V4
    {{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}}, // V5
    {{0, -1, 0}, {1, 0, 0},  {0, 0, 1}},  // V6
    {{0, 1, 0},  {-1, 0, 0}, {0, 0, 1}},  // V7
    {{0, 1, 0},  {1, 0, 0},  {0, 0, -1}}, // V8
    {{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}, // V9
    {{-1, 0, 0}, {0, 0, 1},  {0, 1, 0}},  // V10
    {{1, 0, 0},  {0, 0, -1}, {0, 1, 0}},  // V11
    {{1, 0, 0},  {0, 0, 1},  {0, -1, 0}}, // V12
    {{0, 1, 0},  {0, 0, 1},  {1, 0, 0}},  // V13
    {{0, 1, 0},  {0, 0, -1}, {-1, 0, 0}}, // V14
    {{0, -1, 0}, {0, 0, 1},  {-1, 0, 0}}, // V15
    {{0, -1, 0}, {0, 0, -1}, {1, 0, 0}},  // V16
    {{0, 0, 1},  {1, 0, 0},  {0, 1, 0}},  // V17
    {{0, 0, 1},  {-1, 0, 0}, {0, -1, 0}}, // V18
    {{0, 0, -1}, {1, 0, 0},  {0, -1, 0}}, // V19
    {{0, 0, -1}, {-1, 0, 0}, {0, 1, 0}},  // V20
    {{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}, // V21
    {{0, 0, -1}, {0, 1, 0},  {1, 0, 0}},  // V22
    {{0, 0, 1},  {0, -1, 0}, {1, 0, 0}},  // V23
    {{0, 0, 1},  {0, 1, 0},  {-1, 0, 0}}  // V24
};

int main() {
    // 入力行列を定義
    Matrix3x3 inputMatrixA(3, vector<double>(3));
    cout << endl;
    cout << "一つ目の3x3行列を入力してください (例: 1 0 0 0 1 0 0 0 1 (double型)):" << endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!(cin >> inputMatrixA[i][j])) {
                cerr << "数字以外の入力が検出されました。プログラムを終了します。" << endl;
                return 1;  // 非ゼロの値を返してプログラムを終了
            }
        }
    }

    cout << endl;
    cout << "入力行列は:" << endl;
    printMatrix(inputMatrixA);
    cout << endl;

    try {
        Matrix3x3 inv = inverse(inputMatrixA);
        cout << "逆行列は:" << endl;
        printMatrix(inv);
        cout << endl;
    }
    catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }

    Matrix3x3 inputMatrixB(3, vector<double>(3));
    cout << endl;
    cout << "二つ目の3x3行列を入力してください (例: 1 0 0 0 1 0 0 0 1 (double型)):" << endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!(cin >> inputMatrixB[i][j])) {
                cerr << "数字以外の入力が検出されました。プログラムを終了します。" << endl;
                return 1;  // 非ゼロの値を返してプログラムを終了
            }
        }
    }

    cout << endl;
    cout << "入力行列は:" << endl;
    printMatrix(inputMatrixB);
    cout << endl;

    try {
        Matrix3x3 inv = inverse(inputMatrixB);
        cout << "逆行列は:" << endl;
        printMatrix(inv);
        cout << endl;
    }
    catch (const runtime_error& e) {
        cerr << e.what() << endl;
    }

    // 結果行列を初期化
    Matrix3x3 resultMatrixR(3, vector<double>(3));
    Matrix3x3 invMatrixB = inverse(inputMatrixB);
    multiplyMatrix(inputMatrixA, invMatrixB, resultMatrixR);
    cout << "変換行列Rは:" << endl;
    printMatrix(resultMatrixR);
    cout << endl;

    //// 固有ベクトルを算出
    //cout << "Rの固有ベクトル[H, K, L]は:" << endl;
    //vector<double> b = { 1, 1, 1 };
    //vector<double> x0 = { 0, 0, 0 }; // 初期推定値
    //jacobi(resultMatrixR, b, x0);
    //printVector(x0);
    //cout << endl;


    //// 各固定行列との行列の掛け算と結果の出力
    //for (int index = 0; index < fixedMatrices.size(); ++index) {
    //    multiplyMatrix(inputMatrixA, fixedMatrices[index], resultMatrix);
    //    cout << "固定行列 V" << index + 1 << " との結果行列は:" << endl;
    //    printMatrix(resultMatrix);
    //    cout << endl;
    //}

    return 0;
}
