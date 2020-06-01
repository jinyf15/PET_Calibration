#include "headFile.h"


double pixelprob_det2ML_sum_ph_insert2_double11(double *buf, double *dpos, double *subdpos, double *spos, double *map, int *N_step, double *cos_D, int DM_DET, double mu_det, int fn )
{

	double a, d, sc, x, y,z, pr,  ddx, ddy, dxi, dyi,ddz;
	double sx, sy, sz;	/*/// source position/*/
	double dx, dy, dz;   /*/// detector pixel position/*/
	double px, py, pz;	/*/// projection position on pinhole surface/*/
	double nx,ny,nz,dsq,dsr,dsz,tan_theta;		/*/// incident angle of photon on pinhole surface/*/
	double cx, cy;		/*/// projection of center of ph in the detector surface/*/
	double DDXph, DDYph,DDZ;	/*/// shift displacement in subpixel in ph plane coordinates/*/
	double h, b;
	double theta, beta_phi;	/*///theta cos(beta-phi)/*/
	int i, j, k, l,ii,jj,dpi, kk, ll, indexi, indexj, indexl, indexk, iii, jjj, kkk, lll, NN1, NN2, NN3,NN4;
	double dirt_L, cr,cr3, r2,r3, rho, f,  Rlimit, alpha,t1,i1, tmpi, tmpj, tmpk, tmpl, m1, m2, m3, m4, m, mm, ff, ui[2], uj[2], uk[2], ul[2];
	double prl;
	double w;/*///, *wa;/*/
	double rx,ry,r;
	int Lindex;
	int N1,N2,N3,N4,M1,M2;
	double dirt[4];
	N1 = (int)N_step[0];
	N2 = (int)N_step[1];
	N3 = (int)N_step[2];
	N4 = (int)N_step[3];/*/// change by xiaochun/*/
	NN1 = N1-1;
	NN2 = N2-1;
	NN3 = N3-1;
	NN4 = N4-1;
	M1 = N2*N3*N4;
	M2 = N3*N4;
	dirt[0] = map[N1*N2*N3*N4];
	dirt[1] = map[N1*N2*N3*N4+1];
	dirt[2] = map[N1*N2*N3*N4+2];
	dirt[3] = map[N1*N2*N3*N4+3];/*/// change by xiaochun, the map step/*/

	prl = 0.;
	r2 = buf[5];
	r3 = buf[15];
	Rlimit = buf[9];
	theta = buf[6];
	alpha = buf[10];
	h = buf[11];	/*/// distance between the source and ph surface/*/
	ddx = buf[7];
	ddy = buf[8];
	ddz = buf[13];
	a = fabs(ddx*ddy);

	sx = spos[0];
	sy = spos[1];
	sz = spos[2];
	Lindex = 0;

	for(dpi=0; dpi<DM_DET; dpi++)
	{
		dx = dpos[0+3*dpi];
		dy = dpos[1+3*dpi];
		dz = dpos[2+3*dpi];

		dsr = (dx-sx)*(dx-sx) + (dy-sy)*(dy-sy);
		dsz = (dz-sz)*(dz-sz);
		tan_theta = sqrt(dsr/dsz);
		theta = atan(tan_theta);
		/** pre-exclude if ray cannot pass through the pinhole **/
		if ((theta >= N3*dirt[2]-2*dirt[2]))
		{
			continue;
		}

		ry = -sz/(dz-sz)*(dy-sy)+sy;
	    rx = -sz/(dz-sz)*(dx-sx)+sx;
		r = sqrt(ry*ry + rx*rx);

		/** pre-exclude if outside the rho limit **/
		if (r >= Rlimit){
			continue;
		}
		w = (1-exp(-mu_det*ddz/cos_D[dpi]))*exp((-dpi)*mu_det*ddz/cos_D[dpi]);


		if (1) {
			pr=0.;
			for(ii=0; ii<fn; ii++) {
				for(jj=0; jj<fn; jj++) {
				   // printf("ii %d jj %d\n",ii,jj);
					x = subdpos[Lindex*3+0];
				    y = subdpos[Lindex*3+1];
					z = subdpos[Lindex*3+2];

					ry = -sz/(z-sz)*(y-sy)+sy;
					rx = -sz/(z-sz)*(x-sx)+sx;
					r = sqrt(ry*ry+rx*rx);

                    dsr = (x-sx)*(x-sx) + (y-sy)*(y-sy);
                    dsz = (z-sz)*(z-sz);
                    tan_theta = sqrt(dsr/dsz);
                    theta = atan(tan_theta);
                    tmpk = theta/dirt[2];
                    k = (int)(floor(tmpk));

                    dsq = sqrt(dsr+dsz);
                    nx = (x-sx)/dsq;
                    ny = (y-sy)/dsq;
                    nz = sqrt(nx*nx+ny*ny);
                    nx /= nz;
                    ny /= nz;
                    beta_phi = acos(nx);  // this  needs to be very careful. phi is based on x asix;it has to be consistent with penetration map generation sequence
                    if(ny < 0) beta_phi = 2*PI-beta_phi;

                    tmpl = beta_phi/dirt[3];
                    l = (int)(floor(tmpl));

                    if (l==N4)
                    {
                        tmpl = 0;
                        l = 0;
                    } /*/// ny is nearly equal to 1/*/

					Lindex++;

					if ((r>=Rlimit) || (fabs(rx)>=(N1/2-1)*dirt[0]) || (fabs(ry)>=(N2/2-1)*dirt[1]) || theta >= N3*dirt[2]-2*dirt[2]) { /*/// outside of the rho limit/*/
						continue;
					}
					else{
						tmpi = rx/dirt[0]+N1/2;
						tmpj = ry/dirt[1]+N2/2;
						i = (int)(floor(tmpi));
						j = (int)(floor(tmpj));

						if ((i>=NN1)||(j>=NN2)||(k>=NN3)||(l>=N4)||(i<0)||(j<0)||(k<0)||(l<0)){
							printf("\n <pixelprob_det2ML.c> i %d j %d k %d l %d\n",i,j,k,l);
							printf("N1=%d N2=%d N3=%d N4=%d\n", N1, N2, N3, N4);
							printf("dirt0=%f 1=%f 2=%f 3=%f\n", dirt[0], dirt[1], dirt[2], dirt[3]);
							printf("\n theta %e phi %e",theta,beta_phi);
							printf("\n rx %e ry %e nx %e ny %e ", rx,ry,nx,ny);
							printf("\n r %e Rlimit %e",r,Rlimit);
							printf("\n sx %e %e %e dx %e %e %e",sx,sy,sz,dx,dy,dz);
							printf("\n d %e sc %e",d,sc);
							printf("\n b6 %e b10 %e",buf[6],buf[10]);
						}

						uk[0] = k+1-tmpk;uk[1] = 1-uk[0];
						ul[0] = l+1-tmpl;ul[1] = 1-ul[0];
						ui[0] = i+1-tmpi;ui[1] = 1-ui[0];
						uj[0] = j+1-tmpj;uj[1] = 1-uj[0];
						m = 0;
						for (iii=0; iii<=1; iii++) {
							indexi = i+iii;
							for (jjj=0; jjj<=1; jjj++) {
								indexj = j+jjj;
								for (kkk=0; kkk<=1; kkk++) {
									indexk = k+kkk;
									for (lll=0; lll<=1; lll++) {
										if (l==NN4 && lll==1) indexl = 0;
										else indexl = l+lll;
										m += (map[indexi*M1+indexj*M2+indexk*N4+indexl])*ui[iii]*uj[jjj]*uk[kkk]*ul[lll];
									}
								}
							}
						}

						if (m<0) {
							printf("%f %f %f %f %f %f %f %f\n", uk[0], uk[1], ul[0], ul[1],ui[0], ui[1],uj[0], uj[1]);getchar();
						}
						pr += a*exp(-m);
					}
				}
			}
			prl += pr*cos_D[dpi]/(4*dsq*dsq*PI*fn*fn)*w;
		}
	}		/*/// end of dpi loop/*/
	return (prl);
	/** calculating the prob for partially inculded pixels (end) **/

}
