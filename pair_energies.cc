/*
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <vector>
#include <libmints/mints.h>
#include <libtrans/integraltransform.h>
#define ID(x) ints.DPD_ID(x)

namespace psi { namespace lmp2 {

/* pair_energies(): For RHF references, compute pair energies. Spin-adapt
** pair energies if SPINADAPT_ENERGIES is set to true.
**
** E(IJ) = T2(IJ,AB) * (<ij|ab> - <ij|ba>)
** E(Ij) = T2(Ij,Ab) * <ij|ab>
**
*/

void pair_energies(IntegralTransform& ints, const int& nirrep, const Dimension& occOrbsPI, double** epair)
{
    dpdbuf4 tau, D, E;

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int nab;
    nocc_act = occOrbsPI.sum();
    nab = nocc_act * nocc_act;

    double* pair = new double[nab];

    for(int p=0; p<nab; p++){
        pair[p] = 0.0;
    }

    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, ID("[O>=O]+"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tau, PSIF_CC_TAMPS, 0, ID("[O>=O]+"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "LMP2 tIjAb");
    dpd_buf4_init(&E, PSIF_CC_TAMPS, 0, ID("[O>=O]+"), ID("[O>=O]+"), ID("[O>=O]+"), ID("[O>=O]+"), 0, "E <ij|kl>");
    dpd_contract444(&D, &tau, &E, 0, 0, 1.0, 0.0);
    //dpd_buf4_print(&E, outfile, 1);

    /* Extract diagonal elements (i.e. pair energies) and print them out nicely */
    for(irrep=0; irrep<nirrep; irrep++) {
        double **block;
        dpdparams4 *Params = E.params;
        int p;
        int np = Params->rowtot[irrep];

        dpd_buf4_mat_irrep_init(&E, irrep);
        dpd_buf4_mat_irrep_rd(&E, irrep);
        block = E.matrix[irrep];

        for(p=0; p<np; p++) {
            int i, j, ij;

            i = Params->roworb[irrep][p][0];
            j = Params->roworb[irrep][p][1];

            pair[p] = block[p][p];
        }
        dpd_buf4_mat_irrep_close(&E, irrep);
    }

    *epair = pair;

    dpd_buf4_close(&tau);
    dpd_buf4_close(&D);
    dpd_buf4_close(&E);
}

void print_pair_energies(const int& nirrep, const Dimension& occOrbsPI, double* emp2)
{

    int i, j, ij;
    int irrep;
    int nocc_act = 0;
    int naa, nab;
    nocc_act = occOrbsPI.sum();
    naa = nocc_act * (nocc_act-1)/2;
    nab = nocc_act * nocc_act;

    double emp2_tot = 0.0;

    fprintf(outfile, "\tOrbital pair energies\n");
    fprintf(outfile, "\t    i       j         LMP2\n");
    fprintf(outfile, "\t  -----   -----   ------------\n");
    ij = 0;
    for(i=0; i<nocc_act; i++)
        for(j=0; j<=i; j++,ij++) {
            fprintf(outfile, "\t  %3d     %3d     %12.9lf\n", i+1, j+1, emp2[ij]);
            emp2_tot += emp2[ij];
            if (i != j)
                emp2_tot += emp2[ij];
        }
    fprintf(outfile, "\t  -------------   ------------\n");
    fprintf(outfile, "\t      Total       %12.9lf\n\n", emp2_tot);


    fprintf(outfile, "\n");
}
}} // namespace psi::ccenergy
