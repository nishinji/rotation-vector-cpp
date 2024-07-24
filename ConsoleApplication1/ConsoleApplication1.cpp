#include <iostream>
#include <vector>

using namespace std;

void multiplyMatrix(const vector<vector<float>>& matrixA, const vector<vector<int>>& matrixB, vector<vector<float>>& result) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            result[i][j] = 0.0f;
            for (int k = 0; k < 3; ++k) {
                result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }
}

void printMatrix(const vector<vector<float>>& matrix) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    // 入力行列を定義
    vector<vector<float>> inputMatrix(3, vector<float>(3));
    cout << "3x3行列を入力してください (例: 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0):" << endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            cin >> inputMatrix[i][j];
        }
    }

    // 24個の固定行列を定義
    vector<vector<vector<int>>> fixedMatrices = {
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

    // 結果行列を初期化
    vector<vector<float>> resultMatrix(3, vector<float>(3));

    // 各固定行列との行列の掛け算と結果の出力
    for (int index = 0; index < fixedMatrices.size(); ++index) {
        multiplyMatrix(inputMatrix, fixedMatrices[index], resultMatrix);
        cout << "固定行列 V" << index + 1 << " との結果行列は:" << endl;
        printMatrix(resultMatrix);
        cout << endl;
    }

    return 0;
}
