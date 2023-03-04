//Cost related parameters
float C_day_H = ...; //Daily amortization cost of a HAP 
float C_day_F = ...; //Daily amortization cost of an FSO transceiver on HAP
float C_mtn = ...; //Cost of one time maintenance of a HAP
float D_m = ...; //Maintenance cycle (days)

//Energy parameters
float E_solar = ...; //W-hr
float rho_avion = ...; //W/kg
float rho_HCM_F = ...; //W
float rho_PAT = ...; //W
float rho_inter_F = ...; //W

//HAP-ground link parameters and variables
float sigma = ...;; //3.5*10^-6 (unit 1/m)
float rho_FSO_tx = ...; //W
float R_tx = ...; //m
float rho_rx = ...; //W

float P_tx = ...; //W
float R_rx = ...; //m, Receiver aperture radius

//Other parameters
float H = ...; //m
float L_HH = ...; //m
float W = ...; //The number of wavelengths in WDM technique
float m_H = ...; //kg
float m_F = ...; //kg

//Model parameters
int nHAP = ...; //So luong HAP co san
int nd = ...; //So luong FSO mat dat
int V = ...; //Maximum number of FSO transceivers for interHAP communication
float MAX_XY = ...; //Toa do max cua FSO
int max_M = ...; //So luong max supplementary serving FSO moi HAP
int P = 8; //So canh da giac noi tiep

float a[1..2*P];
float b[1..2*P];
float bigM;
execute {
  for (var k = 1; k <= 2*P; ++k){
    a[k] = 1/Math.tan((k - 0.5) * Math.PI /P);
    b[k] = -Math.cos(Math.PI/(2*P))/Math.sin((k - 0.5) * Math.PI /P);
  }
  bigM = 1/Math.tan(0.5 * Math.PI /P) * MAX_XY + MAX_XY;
}

range rangeHAP = 1..nHAP;
range rangeFSO = 1..nd;
range rangeM = 0..max_M;

float x_FSO[rangeFSO] = ...; //Toa do FSO mat dat
float y_FSO[rangeFSO] = ...; //Toa do FSO mat dat
float R_ext_m[0..max_M] = ...; //Quan he R_ext voi m

dvar float+ x_HAP[rangeHAP]; //Toa do HAP
dvar float+ y_HAP[rangeHAP]; //Toa do HAP
dvar boolean use_m[rangeHAP][rangeM]; //Co su dung m khong
dvar float+ R[rangeHAP];
dvar float useHAP[rangeHAP]; //Co su dung HAPi hay khong
dvar boolean serving[rangeHAP][rangeFSO]; //HAP phuc vu FSO

dexpr float K_over = sum(i in rangeHAP) useHAP[i]; //So luong HAP su dung

//Ham muc tieu
dexpr float Cost_over = K_over * (C_day_H + (1+V) * C_day_F + C_mtn / D_m) + (sum(i in rangeHAP) (sum(j in rangeM) (use_m[i][j] * j))) * C_day_F;
minimize (Cost_over);

subject to
{ 	

  	//Rang buoc giua bien useHAP va serving
 	forall(i in rangeHAP, s in rangeFSO)
 	  	useHAP[i] >= serving[i][s];
 	  	
  	//Moi FSO mat dat phai thuoc mot HAP
  	forall(s in rangeFSO) 
  		sum(i in rangeHAP)(serving[i][s]) == 1;
  		
  	//Moi HAP chi chua toi da W nut mat dat
  	forall(i in rangeHAP)
  	  	sum(s in rangeFSO)(serving[i][s]) <= W;
  	  	
  	//Vi tri nut mat dat nam trong HAP
  	forall(k in 1..P, i in rangeHAP, s in rangeFSO)
  	  serving[i][s] <= 1 - (a[k]*(x_FSO[s] - x_HAP[i]) + (y_FSO[s] - y_HAP[i]) + R[i]*b[k])/bigM;
  	forall(k in P+1..2*P, i in rangeHAP, s in rangeFSO)
  	  serving[i][s] <= 1 - (-a[k]*(x_FSO[s] - x_HAP[i]) - (y_FSO[s] - y_HAP[i]) - R[i]*b[k])/bigM;
 	
 	//Rang buoc su dung m
 	forall(i in rangeHAP)
 	  sum(j in rangeM) use_m[i][j] == 1;
 	  
 	forall(i in rangeHAP)
 	  R[i] <= sum(j in rangeM) use_m[i][j] * R_ext_m[j];
 	  
 	//Cac rang buoc giam khong gian loi giai
 	forall(i in rangeHAP) x_HAP[i] <= MAX_XY;
  	forall(i in rangeHAP) y_HAP[i] <= MAX_XY;
  	
 	forall(i in 2..nHAP) useHAP[i] <= useHAP[i-1];
	forall(i in 1..nHAP - 1) x_HAP[i] >= x_HAP[i + 1];
}