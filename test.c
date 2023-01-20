#define __STDC_WANT_IEC_60559_TYPES_EXT__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#define Jmax 16
#define Kmax 16
#define Lmax 16
#define Mmax 16
typedef _Float128 (*f128_func)(_Float128 *, _Float128 *, _Float128 *);
extern void dqag_sk3d_(_Float128 (*)(_Float128 *, _Float128 *, _Float128 *), _Float128 *xa, _Float128 *xb, _Float128 *ya, _Float128 *yb, _Float128 *za, _Float128 *zb, _Float128 *eps, _Float128 *s, int *ier, int *key);
// extern void dqag_sk3d_(_Float128 (*)(_Float128 *, _Float128 *, _Float128 *), _Float128 *xa, _Float128 *xb, _Float128 *ya, _Float128 *yb, _Float128 *za, _Float128 *zb, _Float128 *eps, _Float128 *s, int *info);
//  _Float128->_Float128にしておいてサジェスト機能を使えるようにする
int main()
{
    int ier = 0;
    int key = 3;
    int info = 0;
    _Float128 xa;
    _Float128 xb;
    _Float128 ya;
    _Float128 yb;
    _Float128 za;
    _Float128 zb;

    char buf[1024];

    //==========物理的な仮定の変数=========
    _Float128 Cap_zeta2 = 100; // 無限遠のzeta
    _Float128 Cap_zetarho = 100;
    // //=========格子点に関する変数=========
    // int Kmax = 16; // 16分割する
    // int Jmax = 16; // 16分割する
    // int Lmax = 16; // 16分割する
    // int Mmax = 16; // 16分割する

    //================変数=============
    _Float128 zeta2[Jmax];
    _Float128 zetarho[Kmax];
    static _Float128 Kjklm[Jmax][Kmax][Lmax][Mmax]; // セグメント違反が出るのでstaticで

    // _Float128 ****Kjklm;
    // Kjklm = (_Float128 ****)malloc(sizeof(int ***) * Jmax);
    // for (int e1 = 0; e1 < Jmax; e1 = e1 + 1)
    // {
    //     Kjklm[e1] = (_Float128 ***)malloc(sizeof(int **) * Kmax); // int*型のスペースを要素数（[][☆][]）の分だけ確保する。
    //     for (int e2 = 0; e2 < Kmax; e2 = e2 + 1)
    //     {
    //         Kjklm[e1][e2] = (_Float128 **)malloc(sizeof(int *) * Lmax); // int型のスペースを要素数（[][][☆]）の分だけ確保する。
    //         for (int e3 = 0; e3 < Lmax; e3 = e3 + 1)
    //         {
    //             Kjklm[e1][e2][e3] = (_Float128 *)malloc(sizeof(int) * Mmax);
    //         }
    //     }
    // }

    // _Float128 Kjklm[Jmax][Kmax][Lmax][Mmax];
    // _Float128 Kjklm[][][][];
    //================

    for (int j = 0; j < Jmax; j++)
    {
        zeta2[j] = Cap_zeta2 * j / Jmax;
    }
    for (int k = 0; k < Kmax; k++)
    {
        zetarho[k] = Cap_zetarho * k / Kmax;
    }
    //===========関数=============
    _Float128 F1_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k)
    {

        return (*zeta2_bar) * (*zeta2_bar) + (*zeta2_j) * (*zeta2_j)           //
               + (*zetarho_bar) * (*zetarho_bar) + (*zetarho_k) * (*zetarho_k) //
               - 2 * ((*zeta2_bar) * (*zeta2_j) + (*zetarho_bar) * (*zetarho_k) * cos(*psi));
    }
    _Float128 F2_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k)
    {

        return ((*zeta2_bar) * (*zeta2_bar) + (*zetarho_bar) * (*zetarho_bar))         //
                   * ((*zeta2_j) * (*zeta2_j) + (*zetarho_k) * (*zetarho_k)) -         //
               ((*zeta2_bar) * (*zeta2_j) + (*zetarho_bar) * (*zetarho_k) * cos(*psi)) //
                   * ((*zeta2_bar) * (*zeta2_j) + (*zetarho_bar) * (*zetarho_k) * cos(*psi));
    }
    _Float128 f_for_psi_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * delta)
    {

        return ((*zeta2_bar) - (*alpha)) * ((*zeta2_bar) - (*beta)) * ((*zetarho_bar) - (*gamma)) * ((*zetarho_bar) - (*delta)) / (((*zeta2_l) - (*alpha)) * ((*zeta2_l) - (*beta)) * ((*zetarho_m) - (*gamma)) * ((*zetarho_m) - (*delta)));
    }
    _Float128 L1(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * delta)

    {
        return 1 / sqrt(F1_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k)) * sqrt(2)                                             //
               / M_PI * (*zetarho_bar) * cos(*psi)                                                                                  //                                                                        //
               * exp(-(*zeta2_bar) * (*zeta2_bar) - (*zetarho_bar) * (*zetarho_bar)                                                 //
                     + F2_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k) / F1_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k)) //
               * f_for_psi_(zeta2_bar, zetarho_bar, zeta2_l, zetarho_m, alpha, beta, gamma, delta);
    }
    _Float128 L2(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * delta)

    {
        return 1 / (2 * sqrt(2) * M_PI) * (*zetarho_bar)                                //
               * sqrt(F1_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k)) * cos(*psi) //
               * exp(-(*zeta2_bar) * (*zeta2_bar) - (*zetarho_bar) * (*zetarho_bar))    //
               * f_for_psi_(zeta2_bar, zetarho_bar, zeta2_l, zetarho_m, alpha, beta, gamma, delta);
    }
    _Float128 Omega1(_Float128 * zeta2_bar, _Float128 * Lambda, _Float128 * alpha_zeta, _Float128 * zeta2_j, _Float128 * zetarho_k, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * de)
    {
        // zetarho_bar = sqrt((*Lambda) * (*Lambda) + 2 * (*zetarho_k) * (*Lambda) * cos(*alpha_zeta) + (*zeta_rho) * (*zeta_rho));
        return sqrt(M_PI / 2) * (*Lambda) * ((*zetarho_k) + (*Lambda) * cos(*alpha_zeta))                                                                                                                                //
               / (sqrt((*Lambda) * (*Lambda) + 2 * (*zetarho_k) * (*Lambda) * cos(*alpha_zeta) + (*zetarho_k) * (*zetarho_k)) * sqrt(((*zeta2_j) - (*zeta2_bar)) * ((*zeta2_j) - (*zeta2_bar)) + (*Lambda) * (*Lambda))) //
               * exp((*zeta2_j) * (*zeta2_j) + (*zetarho_k) * (*zetarho_k)                                                                                                                                               //
                     - ((*zetarho_k) * (*Lambda) * cos(*alpha_zeta) - (*zeta2_j) * ((*zeta2_j) - (*zeta2_bar)) * ((*zeta2_j) - (*zeta2_bar)))                                                                            //
                           * ((*zetarho_k) * (*Lambda) * cos(*alpha_zeta) - (*zeta2_j) * ((*zeta2_j) - (*zeta2_bar)) * ((*zeta2_j) - (*zeta2_bar)))                                                                      //
                           / (((*zeta2_j) - (*zeta2_bar)) * ((*zeta2_j) - (*zeta2_bar)) + (*Lambda) * (*Lambda)));
    }

    //=================================================================
    //============forループ内で使う変数の定義============================
    int j,
        k, l, m;
    _Float128 eps = 1e-8;
    _Float128 s = 0;
    int count = 0;

    //================================================================
    // #pragma omp parallel for private(j, k, l, m, zeta2_j, zetarho_k, zeta2_l, zetarho_m, zeta2_a_l, zetarho_c_m, a_l, b_l, c_m, d_m, alpha, beta, gamma, delta, s, ier) // epsをprivateにすると何故かバグる
    // #pragma omp parallel for private(k, l, m) // firstprivate(zeta2_j, zetarho_k, zeta2_l, zetarho_m, zeta2_a_l, zetarho_c_m, alpha, beta, gamma, delta) // epsをprivateにすると何故かバグる
    // #pragma omp parallel for private(k, l, m, count) // firstprivate(zeta2_j, zetarho_k, zeta2_l, zetarho_m, zeta2_a_l, zetarho_c_m, alpha, beta, gamma, delta) // epsをprivateにすると何故かバグる
    for (j = 0; j < Jmax; j++)
    {
        for (k = 0; k < Kmax; k++)
        {
            for (l = 0; l < Lmax; l++)
            {
                for (m = 0; m < Mmax; m++)
                {
#pragma omp atomic
                    count++;
                    printf("number=%d\n", count);
                    _Float128 psi_a = 0, psi_b = 2 * M_PI;
                    int a_l, b_l, c_m, d_m;                                           // 積分区間のインデックス
                    int a_l_p, a_l_m, b_l_p, b_l_m, c_m_p, c_m_m, d_m_p, d_m_m;       // 偶偶用の積分区間のインデックス(p,mはそれぞれplus,minus,(pp,mm,pm,mp))
                    _Float128 alpha, beta, gamma, delta;                              // psi用の変数
                    _Float128 zeta2_l, zetarho_m, zeta2_j, zetarho_k;                 // zeta[l]とか，特異点
                    _Float128 zeta2_a_l, zeta2_b_l, zetarho_d_m, zetarho_c_m;         // 積分範囲lが始まり，mがおわり
                    _Float128 zeta2_a, zeta2_b, zetarho_d, zetarho_c;                 // 積分範囲lが始まり，mがおわり
                    _Float128 zeta2_a_l_m, zeta2_b_l_m, zetarho_d_m_m, zetarho_c_m_m; // 積分範囲lが始まり，mがおわり

                    _Float128 val_kiki_1, val_kiki_2;
                    _Float128 val_gugu_pp_1, val_gugu_pp_2, val_gugu_pm_1, val_gugu_pm_2, val_gugu_mp_1, val_gugu_mp_2, val_gugu_mm_1, val_gugu_mm_2;
                    _Float128 val_kigu_0p_1, val_kigu_0p_2, val_kigu_0m_1, val_kigu_0m_2, val_guki_p0_1, val_guki_p0_2, val_guki_m0_1, val_guki_m0_2;
                    _Float128 val1, val2, val3, val4; // 特異点を避ける積分範囲の値

                    printf("j= %d,k= %d, l= %d, m= %d\n", j, k, l, m);
                    // printf("singurality"); // ここからprintできない　わけわからん なぜかL1を変えると値も出る

                    zeta2_j = zeta2[j], zetarho_k = zetarho[k], zeta2_l = zeta2[l], zetarho_m = zetarho[m];
                    // printf("singurality");

                    if (j == 0 || k == 0 || k == Kmax)
                    {
                        Kjklm[j][k][l][m] = 0.0;
                    }
                    else
                    {
                        // printf("singurality");

                        if (l % 2 == 1 & m % 2 == 1)
                        {

                            a_l = l - 1, b_l = l + 1, c_m = m - 1, d_m = m + 1;
                            // zeta2_a_l = zeta2[a_l];
                            // #pragma omp critical                                                                                                   // 積分範囲の格子の位置
                            alpha = zeta2[l - 1], beta = zeta2[l + 1], gamma = zetarho[m - 1], delta = zetarho[m + 1]; // psiの0にする点
                                                                                                                       //===========奇数奇数ようの関数
                            _Float128 L1_kiki_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                            {
                                return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                            }
                            _Float128 L2_kiki_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                            {
                                return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                            }

                            if (0 <= a_l & a_l <= Lmax & 0 <= b_l & b_l <= Lmax & 0 <= c_m & c_m <= Mmax & 0 <= d_m & d_m <= Mmax)
                            {
                                // printf("singurality");

                                zeta2_a_l = zeta2[l - 1], zeta2_b_l = zeta2[l + 1], zetarho_c_m = zetarho[m - 1], zetarho_d_m = zetarho[m + 1];

                                if (zeta2[a_l] <= zeta2_j & zeta2_j <= zeta2[b_l] & zetarho[c_m] <= zetarho_k & zetarho_k <= zetarho[d_m])
                                {
                                    // 特異点があるとき
                                    printf("singurality\n");

                                    dqag_sk3d_(&L1_kiki_, &zeta2_a_l, &zeta2_j, &zetarho_c_m, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //                 //左下
                                    val1 = s;
                                    s = 0.0;
                                    // // 右下
                                    dqag_sk3d_(&L1_kiki_, &zeta2_j, &zeta2_b_l, &zetarho_c_m, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //                 //左下

                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kiki_, &zeta2_a_l, &zeta2_j, &zetarho_k, &zetarho_d_m, &psi_a, &psi_b, &eps, &s, &ier, &key); //                 //左下
                                    //  // 左上
                                    val3 = s;
                                    s = 0;
                                    dqag_sk3d_(&L1_kiki_, &zeta2_j, &zeta2_b_l, &zetarho_k, &zetarho_d_m, &psi_a, &psi_b, &eps, &s, &ier, &key); //                 //左下
                                    //  // 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_kiki_1 = val1 + val2 + val3 + val4;
                                    strfromf128(buf, sizeof(buf), "%.40g", val_kiki_1);
                                    printf("val_kiki_2=%s\n", buf);
                                }

                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_kiki_, &zeta2_a_l, &zeta2_b_l, &zetarho_c_m, &zetarho_d_m, &psi_a, &psi_b, &eps, &s, &ier, &key); //                 //左下

                                    val_kiki_1 = s;
                                }
                                printf("L2\n");

                                dqag_sk3d_(&L2_kiki_, &zeta2_a_l, &zeta2_b_l, &zetarho_c_m, &zetarho_d_m, &psi_a, &psi_b, &eps, &s, &ier, &key); //                 //左下

                                val_kiki_2 = s;
                                strfromf128(buf, sizeof(buf), "%.40g", val_kiki_2);
                                printf("val_kiki_2=%s\n", buf);
                                s = 0.0;
                            }
                            else
                            {
                                val_kiki_1 = 0.0, val_kiki_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_kiki_1 + val_kiki_2;
                        }
                        //==================偶数偶数========================================
                        else if (l % 2 == 0 & m % 2 == 0)
                        {
                            a_l_p = l, b_l_p = l + 2, c_m_p = m, d_m_p = m + 2; //// こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                            a_l_m = l - 2, b_l_m = l, c_m_m = m - 2, d_m_m = m;
                            //==============pp==================ppの領域が範囲外に飛び出してなかったら                                                                   // //こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                            if (0 <= a_l_p & a_l_p <= Lmax & 0 <= b_l_p & b_l_p <= Lmax & 0 <= c_m_p & c_m_p <= Mmax & 0 <= d_m_p & d_m_p <= Mmax) // pp
                            {
                                alpha = zeta2[l + 1], beta = zeta2[l + 2], gamma = zetarho[m + 1], delta = zetarho[m + 2];
                                // 関数定義
                                _Float128 L1_gugu_pp(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_gugu_pp(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                } // 間接参照なので並列化によくないかも

                                //===================特異点===============================
                                if (zeta2[a_l_p] <= zeta2_j & zeta2_j <= zeta2[b_l_p] & zetarho[c_m_p] <= zetarho_k & zetarho_k <= zetarho[d_m_p])
                                { // 特異点があるとき;
                                    // 積分範囲
                                    zeta2_a = zeta2[a_l_p], zeta2_b = zeta2[b_l_p], zetarho_c = zetarho[c_m_p], zetarho_d = zetarho[d_m_p];

                                    printf("singurality\n");
                                    dqag_sk3d_(&L1_gugu_pp, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_pp, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_pp, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_pp, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_gugu_pp_1 = val1 + val2 + val3 + val4;
                                    // singurality();
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_gugu_pp, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_gugu_pp_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_gugu_pp, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_gugu_pp_2 = s;
                                s = 0.0;
                            }
                            //==============pp==================ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_pp\n");
                                val_gugu_pp_1 = 0.0;
                                val_gugu_pp_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_gugu_pp_1 + val_gugu_pp_2;
                            //=================pm===================pmの領域が範囲内だったら
                            if (0 <= a_l_p & a_l_p <= Lmax & 0 <= b_l_p & b_l_p <= Lmax & 0 <= c_m_m & c_m_m <= Mmax & 0 <= d_m_m & d_m_m <= Mmax)
                            {
                                alpha = zeta2[l + 1], beta = zeta2[l + 2], gamma = zetarho[m - 1], delta = zetarho[m - 2];
                                _Float128 L1_gugu_pm(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_gugu_pm(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                //=================特異点=======================
                                if (zeta2[a_l_p] <= zeta2_j & zeta2_j <= zeta2[b_l_p] & zetarho[c_m_m] <= zetarho_k & zetarho_k <= zetarho[d_m_m])
                                {
                                    zeta2_a = zeta2[a_l_p], zeta2_b = zeta2[b_l_p], zetarho_c = zetarho[c_m_m], zetarho_d = zetarho[d_m_m];

                                    printf("singurality\n");
                                    dqag_sk3d_(&L1_gugu_pm, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_pm, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_pm, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_pm, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_gugu_pm_1 = val1 + val2 + val3 + val4;
                                    // singurality();
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_gugu_pm, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_gugu_pm_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_gugu_pm, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_gugu_pm_2 = s;
                                s = 0.0;
                            } // ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_pm\n");
                                val_gugu_pm_1 = 0.0;
                                val_gugu_pm_2 = 0.0;
                            }
                            //!-- -- -- -- -- -- -- -- -- -- -- -mp-- -- -- -- -- -- -- -- -- -- -l_m, m_pになる要チェック
                            //
                            if (0 <= a_l_m & a_l_m <= Lmax & 0 <= b_l_m & b_l_m <= Lmax & 0 <= c_m_p & c_m_p <= Mmax & 0 <= d_m_p & d_m_p <= Mmax) // mp
                            {
                                alpha = zeta2[l - 1], beta = zeta2[l - 2], gamma = zetarho[m + 1], delta = zetarho[m + 2];
                                _Float128 L1_gugu_mp(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_gugu_mp(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                if (zeta2[a_l_m] <= zeta2_j & zeta2_j <= zeta2[b_l_m] & zetarho[c_m_p] <= zetarho_k & zetarho_k <= zetarho[d_m_p])
                                {
                                    zeta2_a = zeta2[a_l_m], zeta2_b = zeta2[b_l_m], zetarho_c = zetarho[c_m_p], zetarho_d = zetarho[d_m_p];

                                    printf("singurality\n");
                                    dqag_sk3d_(&L1_gugu_mp, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_mp, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_mp, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_mp, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_gugu_mp_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_gugu_mp, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_gugu_mp_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_gugu_mp, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_gugu_mp_2 = s;
                                s = 0.0;
                            } // ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_mp\n");
                                val_gugu_mp_1 = 0.0;
                                val_gugu_mp_2 = 0.0;
                            }
                            //!-- -- -- -- -- -- -- -- -- -- -- -mm-- -- -- -- -- -- -- -- -- -- -l_m, m_mになる要チェック
                            //
                            if (0 <= a_l_m & a_l_m <= Lmax & 0 <= b_l_m & b_l_m <= Lmax & 0 <= c_m_m & c_m_m <= Mmax & 0 <= d_m_m & d_m_m <= Mmax) // mp
                            {
                                alpha = zeta2[l - 1], beta = zeta2[l - 2], gamma = zetarho[m - 1], delta = zetarho[m - 2];
                                _Float128 L1_gugu_mm(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_gugu_mm(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                if (zeta2[a_l_m] <= zeta2_j & zeta2_j <= zeta2[b_l_m] & zetarho[c_m_m] <= zetarho_k & zetarho_k <= zetarho[d_m_m])
                                {
                                    zeta2_a = zeta2[a_l_m], zeta2_b = zeta2[b_l_m], zetarho_c = zetarho[c_m_m], zetarho_d = zetarho[d_m_m];

                                    printf("singurality\n");
                                    dqag_sk3d_(&L1_gugu_mm, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_mm, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_mm, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_gugu_mm, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_gugu_mm_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_gugu_mm, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_gugu_mm_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_gugu_mm, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_gugu_mm_2 = s;
                                s = 0.0;
                            } // ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_mm\n");
                                val_gugu_mm_1 = 0.0;
                                val_gugu_mm_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_gugu_pp_1 + val_gugu_pp_2 + val_gugu_pm_1 + val_gugu_pm_2 //
                                                + val_gugu_mp_1 + val_gugu_mp_2 + val_gugu_mm_1 + val_gugu_mm_2;
                        }

                        //==================奇数偶数====================
                        else if (l % 2 == 1 & m % 2 == 0)
                        {
                            a_l = l - 1, b_l = l + 1, c_m_p = m, d_m_p = m + 2, c_m_m = m - 2, d_m_m = m;
                            //=============0p=========================
                            if (0 <= a_l & a_l <= Lmax & 0 <= b_l & b_l <= Lmax & 0 <= c_m_p & c_m_p <= Mmax & 0 <= d_m_p & d_m_p <= Mmax)
                            {
                                alpha = zeta2[l - 1], beta = zeta2[l + 1], gamma = zetarho[m + 1], delta = zetarho[m + 2];
                                _Float128 L1_kigu_0p(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_kigu_0p(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                if (zeta2[a_l] <= zeta2_j & zeta2_j <= zeta2[b_l] & zetarho[c_m_p] <= zetarho_k & zetarho_k <= zetarho[d_m_p])
                                {
                                    printf("singrality");
                                    dqag_sk3d_(&L1_kigu_0p, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kigu_0p, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kigu_0p, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kigu_0p, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_kigu_0p_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_kigu_0p, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_kigu_0p_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_kigu_0p, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_kigu_0p_2 = s;
                                s = 0.0;
                            }
                            else
                            {
                                printf("hanigai_0p\n");
                                val_kigu_0p_1 = 0.0;
                                val_kigu_0p_2 = 0.0;
                            }

                            //================0m=======================
                            if (0 <= a_l & a_l <= Lmax & 0 <= b_l & b_l <= Lmax & 0 <= c_m_m & c_m_m <= Mmax & 0 <= d_m_m & d_m_m <= Mmax)
                            {
                                alpha = zeta2[l - 1], beta = zeta2[l + 1], gamma = zetarho[m - 1], delta = zetarho[m - 2];
                                _Float128 L1_kigu_0m(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_kigu_0m(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                if (zeta2[a_l] <= zeta2_j & zeta2_j <= zeta2[b_l] & zetarho[c_m_p] <= zetarho_k & zetarho_k <= zetarho[d_m_p])
                                {
                                    printf("singrality");
                                    dqag_sk3d_(&L1_kigu_0m, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kigu_0m, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kigu_0m, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_kigu_0m, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_kigu_0m_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_kigu_0m, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_kigu_0m_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_kigu_0m, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_kigu_0m_2 = s;
                                s = 0.0;
                            }
                            else
                            {
                                printf("hanigai_0m\n");
                                val_kigu_0m_1 = 0.0;
                                val_kigu_0m_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_kigu_0p_1 + val_kigu_0p_2 + val_kigu_0m_1 + val_kigu_0m_2;
                        }

                        //===============偶数奇数===========
                        else
                        {
                            a_l_p = l, b_l_p = l + 2, a_l_m = l - 2, b_l_p = l, c_m = m - 1, d_m = m + 1;
                            //================p0===============
                            if (0 <= a_l_p & a_l_p <= Lmax & 0 <= b_l_p & b_l_p <= Lmax & 0 <= c_m & c_m <= Mmax & 0 <= d_m & d_m <= Mmax)
                            {
                                alpha = zeta2[l + 1], beta = zeta2[l + 2], gamma = zetarho[m + 1], delta = zetarho[m - 1];
                                _Float128 L1_guki_p0(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_guki_p0(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }

                                //======特異点があるとき=========
                                if (zeta2[a_l_p] <= zeta2_j & zeta2_j <= zeta2[b_l_p] & zetarho[c_m] <= zetarho_k & zetarho_k <= zetarho[d_m])
                                {
                                    printf("singrality");
                                    dqag_sk3d_(&L1_guki_p0, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_guki_p0, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_guki_p0, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_guki_p0, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_guki_p0_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_guki_p0, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_guki_p0_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_guki_p0, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_guki_p0_2 = s;
                                s = 0.0;
                            }
                            else
                            {
                                printf("hanigai_p0\n");
                                val_guki_p0_1 = 0.0;
                                val_guki_p0_2 = 0.0;
                            }
                            //================m0================
                            if (0 <= a_l_m & a_l_m <= Lmax & 0 <= b_l_m & b_l_m <= Lmax & 0 <= c_m & c_m <= Mmax & 0 <= d_m & d_m <= Mmax)
                            {
                                alpha = zeta2[l - 1], beta = zeta2[l - 2], gamma = zetarho[m + 1], delta = zetarho[m - 1];
                                _Float128 L1_guki_m0(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L1(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                _Float128 L2_guki_m0(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi)
                                {
                                    return L2(zeta2_bar, zetarho_bar, psi, &zeta2_j, &zetarho_k, &zeta2_l, &zetarho_m, &alpha, &beta, &gamma, &delta);
                                }
                                //===========特異点があるとき===========
                                if (zeta2[a_l_m] <= zeta2_j & zeta2_j <= zeta2[b_l_m] & zetarho[c_m] <= zetarho_k & zetarho_k <= zetarho[d_m])
                                {
                                    printf("singrality");
                                    dqag_sk3d_(&L1_guki_m0, &zeta2_a, &zeta2_j, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); // 左下;
                                    val1 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_guki_m0, &zeta2_j, &zeta2_b, &zetarho_c, &zetarho_k, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右下
                                    val2 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_guki_m0, &zeta2_a, &zeta2_j, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 左上
                                    val3 = s;
                                    s = 0.0;
                                    dqag_sk3d_(&L1_guki_m0, &zeta2_j, &zeta2_b, &zetarho_k, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key); //; 右上
                                    val4 = s;
                                    s = 0.0;
                                    val_guki_m0_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    dqag_sk3d_(&L1_guki_m0, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                    val_guki_m0_1 = s;
                                    s = 0.0;
                                }
                                printf("L2\n");
                                dqag_sk3d_(&L2_guki_m0, &zeta2_a, &zeta2_b, &zetarho_c, &zetarho_d, &psi_a, &psi_b, &eps, &s, &ier, &key);
                                val_guki_m0_2 = s;
                                s = 0.0;
                            }
                            else
                            {
                                printf("hanigai_m0\n");
                                val_guki_m0_1 = 0.0;
                                val_guki_m0_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_guki_p0_1 + val_guki_p0_2 + val_guki_m0_1 + val_guki_m0_2;
                        }
                    }
                    strfromf128(buf, sizeof(buf), "%.40g", Kjklm[j][k][l][m]);
                    printf("s=%s\n", buf);
                }
            }
        }
    }
    FILE *file;
    file = fopen("test.dat", "wb");
    fwrite(Kjklm, sizeof(Kjklm), 1, file);
    fclose(file);
    return 0;
}