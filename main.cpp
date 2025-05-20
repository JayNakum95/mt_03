#include "Novice.h"
#include <cmath>
#include <cstdio>

const char kWindowTitle[] = "GC2B_07_ナクム_ジェイ_ハルシュバルダン";

// 行列・ベクトル定義
struct Matrix4x4 {
    float m[4][4];
};
struct Matrix3x3 {
    float m[3][3];
};
struct Vector3 {
    float x, y, z;
};


// 画面表示用定数
static const int kRowHeight = 20;
static const int kColumnWidth = 60;

// 行列の加算
/*static Matrix4x4 Add(Matrix4x4& m1, Matrix4x4& m2) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            result.m[row][col] = m1.m[row][col] + m2.m[row][col];
    return result;
}

// 行列の減算
static Matrix4x4 Subtract(Matrix4x4& m1, Matrix4x4& m2) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            result.m[row][col] = m1.m[row][col] - m2.m[row][col];
    return result;
}

// 行列の乗算
static Matrix4x4 Multiply(Matrix4x4& m1, Matrix4x4& m2) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            result.m[row][col] = 0;
            for (int k = 0; k < 4; ++k) {
                result.m[row][col] += m1.m[row][k] * m2.m[k][col];
            }
        }
    }
    return result;
}

// 行列の転置
static Matrix4x4 Transpose(const Matrix4x4& m) {
    Matrix4x4 result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            result.m[row][col] = m.m[col][row];
    return result;
}
*/
// 単位行列
static Matrix4x4 MakeIdentity4x4() {
    Matrix4x4 result{};
    for (int i = 0; i < 4; ++i)
        result.m[i][i] = 1.0f;
    return result;
}

// 行列の逆行列（簡易版・完全ではない）
// 4x4行列の逆行列（一般的なガウス・ジョルダン法）
/*/static Matrix4x4 Inverse(Matrix4x4& m) {
    Matrix4x4 result = MakeIdentity4x4();
    Matrix4x4 temp = m;

    for (int i = 0; i < 4; ++i) {
        // 対角成分が0の場合、行を入れ替える
        if (fabsf(temp.m[i][i]) < 1e-6f) {
            bool swapped = false;
            for (int j = i + 1; j < 4; ++j) {
                if (fabsf(temp.m[j][i]) > 1e-6f) {
                    for (int k = 0; k < 4; ++k) {
                        std::swap(temp.m[i][k], temp.m[j][k]);
                        std::swap(result.m[i][k], result.m[j][k]);
                    }
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                // 逆行列が存在しない
                return MakeIdentity4x4(); // 代替措置
            }
        }

        // 対角要素を1にする
        float diag = temp.m[i][i];
        for (int k = 0; k < 4; ++k) {
            temp.m[i][k] /= diag;
            result.m[i][k] /= diag;
        }

        // 他の行のi列を0にする
        for (int j = 0; j < 4; ++j) {
            if (i == j) continue;
            float factor = temp.m[j][i];
            for (int k = 0; k < 4; ++k) {
                temp.m[j][k] -= factor * temp.m[i][k];
                result.m[j][k] -= factor * result.m[i][k];
            }
        }
    }

    return result;
}*/
Matrix4x4 MakeTranslateMatrix(const Vector3& translate) {
    Matrix4x4 result = MakeIdentity4x4();
    result.m[0][3] = translate.x;
    result.m[1][3] = translate.y;
    result.m[2][3] = translate.z;
    return result;
}
Matrix4x4 MakeScaleMatrix(const Vector3& scale) {
    Matrix4x4 result{};
    result.m[0][0] = scale.x;
    result.m[1][1] = scale.y;
    result.m[2][2] = scale.z;
    result.m[3][3] = 1.0f;
    return result;
}
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
    Vector3 result;
    // Define Matrix3x3 structure to resolve the undefined identifier error

    result.x = vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0];
    result.y = vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1];
    result.z = vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2];
    float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

    // Avoid division by zero
    if (w != 0.0f) {
        result.x /= w;
        result.y /= w;
        result.z /= w;
    }

    return result;
}

void VectorScreenPrintf(int x, int y, const Vector3& vector, const char* label) {
    Novice::ScreenPrintf(x, y, "%s: (%.2f, %.2f, %.2f)", label, vector.x, vector.y, vector.z);
}



// 画面に行列を表示
static void MatrixScreenPrintf(int x, int y, Matrix4x4& matrix, const char* label) {
    Novice::ScreenPrintf(x, y - 20, "%s", label);
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            Novice::ScreenPrintf(
                x + column * kColumnWidth,
                y + row * kRowHeight,
                "%6.02f",
                matrix.m[row][column]);
        }
    }
}

// メイン
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
    Novice::Initialize(kWindowTitle, 1280, 720);

    Matrix4x4 m1 = {
        3.2f, 0.7f, 9.6f, 4.4f,
        5.5f, 1.3f, 7.8f, 2.1f,
        6.9f, 8.0f, 2.6f, 1.0f,
        0.5f, 7.2f, 5.1f, 3.3f
    };

    Matrix4x4 m2 = {
        4.1f, 6.5f, 3.3f, 2.2f,
        8.8f, 0.6f, 9.9f, 7.7f,
        1.1f, 5.5f, 6.6f, 0.0f,
        3.3f, 9.9f, 8.8f, 2.2f
    };

    // 計算
   /* Matrix4x4 resultAdd = Add(m1, m2);
    Matrix4x4 resultSubtract = Subtract(m1, m2);
    Matrix4x4 resultMultiply = Multiply(m1, m2);
    Matrix4x4 inverseM1 = Inverse(m1);  // 仮の単位行列
    Matrix4x4 inverseM2 = Inverse(m2);  // 仮の単位行列
    Matrix4x4 transposeM1 = Transpose(m1);
    Matrix4x4 transposeM2 = Transpose(m2);
    Matrix4x4 identity = MakeIdentity4x4();*/
    Vector3 translate{ 4.1f, 2.6f, 0.8f };
    Vector3 scale{ 1.5f, 5.2f, 7.3f };
    Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);
    Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
    Vector3 point{ 2.3f, 3.8f, 1.4f };
    Matrix4x4 transformMatrix = {
        1.0f, 2.0f, 3.0f, 4.0f,
        3.0f, 1.0f, 1.0f, 2.0f,
        1.0f, 4.0f, 2.0f, 3.0f,
        2.0f, 2.0f, 1.0f, 3.0f,

    };
    // ループ
    while (Novice::ProcessMessage() == 0) {
        Novice::BeginFrame();

        /* MatrixScreenPrintf(0, 20, resultAdd, "Add");
         MatrixScreenPrintf(0, 20 + kRowHeight * 5, resultSubtract, "Subtract");
         MatrixScreenPrintf(0, 20 + kRowHeight * 5 * 2, resultMultiply, "Multiply");
         MatrixScreenPrintf(0, 20 + kRowHeight * 5 * 3, inverseM1, "inverseM1");
         MatrixScreenPrintf(0, 20 + kRowHeight * 5 * 4, inverseM2, "inverseM2");
         MatrixScreenPrintf(kColumnWidth * 5, 20, transposeM1, "transposeM1");
         MatrixScreenPrintf(kColumnWidth * 5, 20 + kRowHeight * 5, transposeM2, "transposeM2");
         MatrixScreenPrintf(kColumnWidth * 5, 20 + kRowHeight * 5 * 2, identity, "identity");*/
        Vector3 transformed = Transform(point, transformMatrix);
        VectorScreenPrintf(0, 0, transformed, "transformed");
        MatrixScreenPrintf(0, 40, translateMatrix, "translateMatrix");
        MatrixScreenPrintf(0, 40 + kRowHeight * 5, scaleMatrix, "scaleMatrix");
        Novice::EndFrame();
    }

    Novice::Finalize();
    return 0;
}
