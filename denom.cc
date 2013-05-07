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
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libmints/mints.h>

namespace psi { namespace lmp2 {

void denom(const int& nirrep,
           const Dimension& occOrbsPI,
           const Dimension& virOrbsPI,
           const Dimension& occOffset,
           const Dimension& virOffset)
{
    Dimension openpi(nirrep);

    int h, i, j, a, b, ij, ab;
    int I, J, A, B;
    int isym, jsym, asym, bsym;
    double fii, fjj, faa, fbb;
    dpdfile2 fIJ, fij, fAB, fab;
    dpdfile2 dIA, dia;
    dpdfile4 dIJAB, dijab, dIjAb;

    /* Grab Fock matrices from disk */
    dpd_file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_mat_init(&fIJ);
    dpd_file2_mat_rd(&fIJ);

    dpd_file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    dpd_file2_mat_init(&fij);
    dpd_file2_mat_rd(&fij);

    dpd_file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_mat_init(&fAB);
    dpd_file2_mat_rd(&fAB);

    dpd_file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    dpd_file2_mat_init(&fab);
    dpd_file2_mat_rd(&fab);

    /* Alpha one-electron denominator */
    dpd_file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_mat_init(&dIA);

    for(h=0; h < nirrep; h++) {

        for(i=0; i < occOrbsPI[h]; i++) {
            fii = fIJ.matrix[h][i][i];

            for(a=0; a < (virOrbsPI[h] - openpi[h]); a++) {
                faa = fAB.matrix[h][a][a];

                dIA.matrix[h][i][a] = 1.0/(fii - faa);
            }
        }
    }

    dpd_file2_mat_wrt(&dIA);
    dpd_file2_mat_close(&dIA);
    dpd_file2_close(&dIA);

    /* Beta one-electron denominator */
    dpd_file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
    dpd_file2_mat_init(&dia);

    for(h=0; h < nirrep; h++) {

        for(i=0; i < (occOrbsPI[h] - openpi[h]); i++) {
            fii = fij.matrix[h][i][i];

            for(a=0; a < virOrbsPI[h]; a++) {
                faa = fab.matrix[h][a][a];

                dia.matrix[h][i][a] = 1.0/(fii - faa);
            }
        }
    }

    dpd_file2_mat_wrt(&dia);
    dpd_file2_mat_close(&dia);
    dpd_file2_close(&dia);

    /* Alpha-alpha two-electron denominator */
    dpd_file4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, "dIJAB");

    for(h=0; h < nirrep; h++) {

        dpd_file4_mat_irrep_init(&dIJAB, h);

        /* Loop over the rows */
        for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
            i = dIJAB.params->roworb[h][ij][0];
            j = dIJAB.params->roworb[h][ij][1];
            isym = dIJAB.params->psym[i];
            jsym = dIJAB.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occOffset[isym];
            J = j - occOffset[jsym];

            fii = fIJ.matrix[isym][I][I];
            fjj = fIJ.matrix[jsym][J][J];

            /* Loop over the columns */
            for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
                a = dIJAB.params->colorb[h][ab][0];
                b = dIJAB.params->colorb[h][ab][1];
                asym = dIJAB.params->rsym[a];
                bsym = dIJAB.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - virOffset[asym];
                B = b - virOffset[bsym];

                faa = fAB.matrix[asym][A][A];
                fbb = fAB.matrix[bsym][B][B];

                dIJAB.matrix[h][ij][ab] =
                    ((A >= (virOrbsPI[asym] - openpi[asym])) ||
                     (B >= (virOrbsPI[bsym] - openpi[bsym])) ?
                     0.0 : 1.0/(fii + fjj - faa - fbb));
            }
        }

        dpd_file4_mat_irrep_wrt(&dIJAB, h);
        dpd_file4_mat_irrep_close(&dIJAB, h);

    }

    dpd_file4_close(&dIJAB);

    /* Beta-beta two-electron denominator */
    dpd_file4_init(&dijab, PSIF_CC_DENOM, 0, 1, 6, "dijab");

    for(h=0; h < nirrep; h++) {

        dpd_file4_mat_irrep_init(&dijab, h);

        /* Loop over the rows */
        for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
            i = dijab.params->roworb[h][ij][0];
            j = dijab.params->roworb[h][ij][1];
            isym = dijab.params->psym[i];
            jsym = dijab.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occOffset[isym];
            J = j - occOffset[jsym];

            fii = fij.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            /* Loop over the columns */
            for(ab=0; ab < dijab.params->coltot[h]; ab++) {
                a = dijab.params->colorb[h][ab][0];
                b = dijab.params->colorb[h][ab][1];
                asym = dijab.params->rsym[a];
                bsym = dijab.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - virOffset[asym];
                B = b - virOffset[bsym];

                faa = fab.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dijab.matrix[h][ij][ab] =
                    ((I >= (occOrbsPI[isym] - openpi[isym])) ||
                     (J >= (occOrbsPI[jsym] - openpi[jsym])) ?
                     0.0 : 1.0/(fii + fjj - faa - fbb));
            }
        }

        dpd_file4_mat_irrep_wrt(&dijab, h);
        dpd_file4_mat_irrep_close(&dijab, h);

    }

    dpd_file4_close(&dijab);

    /* Alpha-beta two-electron denominator */
    dpd_file4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, "dIjAb");

    for(h=0; h < nirrep; h++) {

        dpd_file4_mat_irrep_init(&dIjAb, h);

        /* Loop over the rows */
        for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
            i = dIjAb.params->roworb[h][ij][0];
            j = dIjAb.params->roworb[h][ij][1];
            isym = dIjAb.params->psym[i];
            jsym = dIjAb.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occOffset[isym];
            J = j - occOffset[jsym];

            fii = fIJ.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            /* Loop over the columns */
            for(ab=0; ab < dIjAb.params->coltot[h]; ab++) {
                a = dIjAb.params->colorb[h][ab][0];
                b = dIjAb.params->colorb[h][ab][1];
                asym = dIjAb.params->rsym[a];
                bsym = dIjAb.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - virOffset[asym];
                B = b - virOffset[bsym];

                faa = fAB.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dIjAb.matrix[h][ij][ab] =
                    ((A >= (virOrbsPI[asym] - openpi[asym])) ||
                     (J >= (occOrbsPI[jsym] - openpi[jsym])) ?
                     0.0 : 1.0/(fii + fjj - faa - fbb));
            }
        }

        dpd_file4_mat_irrep_wrt(&dIjAb, h);
        dpd_file4_mat_irrep_close(&dIjAb, h);

    }

    dpd_file4_close(&dIjAb);

    dpd_file2_mat_close(&fIJ);
    dpd_file2_mat_close(&fij);
    dpd_file2_mat_close(&fAB);
    dpd_file2_mat_close(&fab);
    dpd_file2_close(&fIJ);
    dpd_file2_close(&fij);
    dpd_file2_close(&fAB);
    dpd_file2_close(&fab);

}

}} // namespace psi::ccsort
