from math import pi, sqrt
from cmath import exp
from matplotlib import pyplot as plt
import numpy.linalg as la
import numpy as np

# t - NN hopping # tL - IL hopping

class blg_site :
	def __init__(self, name, unit_x, unit_y, site_eps_fn) :
		self.name = name
		
		self.unit_x = unit_x
		self.unit_y = unit_y
		
		self.x = unit_x
		if self.name == 'AU' :
			self.y = unit_y + 1/sqrt(3)
		elif self.name == 'BL' :
			self.y = unit_y - 1/sqrt(3)
		else :
			self.y = unit_y 

		self.site_eps = site_eps_fn(self.name, self.x, self.y)

	def __str__(self):
		return self.name + ' ' + str(self.unit_x) + ' ' + str(self.unit_y)

	def onsite(self, param) :
		return self.site_eps + param	
	
	def unit_cell_pos(self) :
		return np.array([self.unit_x,self.unit_y])

	def hopping(self, site ,param) : # param[0] = t # param[1] = tL
		if abs(self.unit_x-site.unit_x) < 1e-7 and abs(self.unit_y-site.unit_y) < 1e-7 : 
			if self.name == 'AU' and site.name == 'BU' :
				return param[0]
			if self.name == 'BU' and site.name == 'AU' :
				return param[0]
			if self.name == 'BU' and site.name == 'AL' :
				return param[1]
			if self.name == 'AL' and site.name == 'BU' :
				return param[1]
			if self.name == 'AL' and site.name == 'BL' :
				return param[0]
			if self.name == 'BL' and site.name == 'AL' :
				return param[0]
		
		if abs(self.unit_x-site.unit_x+1/2) < 1e-7 and abs(self.unit_y-site.unit_y+sqrt(3)/2) < 1e-7 :	
		# if self.unit_x == site.unit_x - 1/2 and self.unit_y == site.unit_y - sqrt(3)/2 : 
			if self.name == 'AU' and site.name == 'BU' :
				return param[0]
			if self.name == 'AL' and site.name == 'BL' :
				return param[0]
		if abs(self.unit_x-site.unit_x+1/2) < 1e-7 and abs(self.unit_y-site.unit_y-sqrt(3)/2) < 1e-7 :
		# if self.unit_x == site.unit_x - 1/2 and self.unit_y == site.unit_y + sqrt(3)/2 :
			if self.name == 'BU' and site.name == 'AU' :
				return param[0]
			if self.name == 'BL' and site.name == 'AL' :
				return param[0]
		if abs(self.unit_x-site.unit_x-1/2) < 1e-7 and abs(self.unit_y-site.unit_y-sqrt(3)/2) < 1e-7 :
		# if self.unit_x == site.unit_x + 1/2 and self.unit_y == site.unit_y + sqrt(3)/2 :
			if self.name == 'BU' and site.name == 'AU' :
				return param[0]
			if self.name == 'BL' and site.name == 'AL' :
				return param[0]
		
		if abs(self.unit_x-site.unit_x-1/2) < 1e-7 and abs(self.unit_y-site.unit_y+sqrt(3)/2) < 1e-7 :
		# if self.unit_x == site.unit_x + 1/2 and self.unit_y == site.unit_y - sqrt(3)/2 :
			if self.name == 'AU' and site.name == 'BU' :
				return param[0]
			if self.name == 'AL' and site.name == 'BL' :
				return param[0]
		
		return 0


class blg_unit_cell :
	def __init__(self, x, y, site_eps_fn) : 
		self.x = x
		self.y = y

		self.sites = []
		self.sites.append(blg_site('BL',x,y,site_eps_fn))
		self.sites.append(blg_site('BU',x,y,site_eps_fn))
		self.sites.append(blg_site('AL',x,y,site_eps_fn))
		self.sites.append(blg_site('AU',x,y,site_eps_fn))
	
	def __str__(self) :
		to_str = str(list(map(str ,self.sites)))
		return to_str

	def remove(self,name) :
		for site in self.sites :
			if site.name == name :
				self.sites.remove(site)



class blg_subcolumn1 :
	def __init__(self, x, W, site_eps_fn) :
		self.num_unit_cells = int(W/sqrt(3))
		
		self.unit_cells = []
		for i in range(self.num_unit_cells) :
			y = i*sqrt(3) + sqrt(3)/2
			self.unit_cells.append(blg_unit_cell(x, y, site_eps_fn))

		self.sites = []
		for i in range(self.num_unit_cells) :
			self.sites = self.sites + self.unit_cells[i].sites

	def hamiltonian(self,U,t,tL) :
		N = len(self.sites)
		H = np.zeros((N,N))
		
		for i in range(N) :
			for j in range(N) :
				if i == j :
					H[i,j] = self.sites[i].onsite(U) 		
				else :
					H[i,j] = self.sites[i].hopping(self.sites[j],[-t,-tL])

		return H	




class blg_subcolumn2 :
	def __init__(self, x, W, site_eps_fn) :
		self.num_unit_cells = int(W/sqrt(3)) + 1
		
		self.unit_cells = []
		for i in range(self.num_unit_cells) :
			y = i*sqrt(3)
			self.unit_cells.append(blg_unit_cell(x, y, site_eps_fn))

		self.unit_cells[0].remove('BU')
		self.unit_cells[0].remove('BL')
		self.unit_cells[-1].remove('AU')
		self.unit_cells[-1].remove('AL')	
		self.sites = []
		for i in range(self.num_unit_cells) :
			self.sites = self.sites + self.unit_cells[i].sites

	def hamiltonian(self,U,t,tL) :
		N = len(self.sites)
		H = np.zeros((N,N))
		
		for i in range(N) :
			for j in range(N) :
				if i == j :
					H[i,j] = self.sites[i].onsite(U) 		
				else :
					H[i,j] = self.sites[i].hopping(self.sites[j],[-t,-tL])

		return H




class blg_lead :
	def __init__(self, W, site_eps_fn,t,tL) :
		
		self.t = t
		self.tL = tL
		self.W = W	
		self.subcolumns = []
		self.subcolumns.append(blg_subcolumn1(0,W,site_eps_fn))
		self.subcolumns.append(blg_subcolumn2(0.5,W,site_eps_fn))
		self.sites = self.subcolumns[0].sites + self.subcolumns[1].sites
		self.debug = 0

	def hamiltonian(self) :	
		N = len(self.sites)
		H = np.zeros((N,N))
		
		for i in range(N) :
			for j in range(N) :
				if i == j :
					H[i,j] = self.sites[i].onsite(0) 		
				else :
					H[i,j] = self.sites[i].hopping(self.sites[j],[-self.t,-self.tL])
		
		n = int(N/2)
		self.debug = H
		H11 = H[:n,:n]
		H22 = H[n:,n:]
		H12 = H[:n,n:]
		H21 = H[n:,:n]

		return H11, H12, H21, H22
		
	def modes(self,E) :
		H11, H12, H21, H22 = self.hamiltonian()
		N = int(len(self.sites)/2)

		M11 = -np.eye(N);
		M12 = la.inv(H21) @ (E*np.eye(N)-H22);
		M21 = -1*la.inv(H12) @ (E*np.eye(N)-H11);
		M22 = (la.inv(H12) @ (E*np.eye(N)-H11) @ la.inv(H21) @ (E*np.eye(N)-H22)) - np.eye(N);		
		M = np.block([[M11, M12],[M21, M22]])

		w,v = la.eig(M)

		return w,v,M

	def vel_fn(self) :
		H11, H12, H21, H22 = self.hamiltonian()
		O = np.zeros(np.shape(H11))
		beta_dag = np.block([[O, H12],[O, O]])
		beta = np.block([[O, O],[H21, O]])	
		return lambda k: 1j*exp(1j*k)*beta - 1j*exp(-1j*k)*beta_dag

	def positive_modes(self,E) :
		w,v,M = self.modes(E)
		kw = np.log(w)/1j
		N = len(w)
		wp = []
		vp = []
		for i in range(N) :
			if abs(np.imag(kw[i])) > 1e-10 and np.imag(kw[i]) > 0 :
				wp.append(w[i])
				vp.append(v[:,i])

			if abs(np.imag(kw[i])) <= 1e-10 :
				vf = self.vel_fn()
				vel = np.conjugate(np.transpose(v[:,i:i+1]))@vf(kw[i])@v[:,i:i+1]
				if vel[0] > 0 :
					wp.append(w[i])
					vp.append(v[:,i])					

		wp = np.array(wp)
		vp = np.transpose(np.array(vp))
		return wp,vp

	def negative_modes(self,E) :
		w,v,M = self.modes(E)
		kw = np.log(w)/1j
		N = len(w)
		wn = []
		vn = []
		for i in range(N) :
			if abs(np.imag(kw[i])) > 1e-10 and np.imag(kw[i]) <= 0 :
				wn.append(w[i])
				vn.append(v[:,i])

			if abs(np.imag(kw[i])) <= 1e-10 :
				vf = self.vel_fn()
				vel = np.conjugate(np.transpose(v[:,i:i+1]))@vf(kw[i])@v[:,i:i+1]
				if vel[0] <= 0 :
					wn.append(w[i])
					vn.append(v[:,i])					

		wn = np.array(wn)
		vn = np.transpose(np.array(vn))
		return wn,vn	

	def plot_bandstructure(self,k_arr) :
		H11, H12, H21, H22 = self.hamiltonian()
		O = np.zeros(np.shape(H11))
		beta_dag = np.block([[O, H12],[O, O]])
		beta = np.block([[O, O],[H21, O]])
		alpha = np.block([[H11, H12],[H21, H22]])
		H_fun = lambda k: alpha + exp(1j*k)*beta + exp(-1j*k)*beta_dag
		Ek = []
		for k in k_arr :
			val = la.eigvals(H_fun(k))
			Ek.append(sorted(np.real(val)))
		Ek = np.array(Ek)
		plt.figure()
		plt.plot(k_arr,Ek)
		plt.show()

	def subcolumn_modes(self,E) :
		wp,vp = self.positive_modes(E)
		wn,vn = self.negative_modes(E)
		N = len(wp)
		Up1 = vp[:N,:]
		Up2 = vp[N:,:]
		Un1 = vn[:N,:]
		Un2 = vn[N:,:]

		return Up1, Up2, Un1, Un2, wp, wn

	def bloch_matrices(self,E) :
		Up1, Up2, Un1, Un2, wp, wn = self.subcolumn_modes(E)
		Fp1 = Up1@np.diag(wp)@la.inv(Up1)
		Fp2 = Up2@np.diag(wp)@la.inv(Up2)
		Fn1 = Un1@np.diag(wn)@la.inv(Un1)
		Fn2 = Un2@np.diag(wn)@la.inv(Un2)

		return Fp1, Fp2, Fn1, Fn2

	def self_energy(self,E) :
		Up1, Up2, Un1, Un2, wp, wn = self.subcolumn_modes(E)
		H11, H12, H21, H22 = self.hamiltonian()
		Fp1, Fp2, Fn1, Fn2 = self.bloch_matrices(E)
		N = len(wp)
		ML = H12@la.inv(E*np.eye(N)-H22)@H21@(la.inv(Fn1)+np.eye(N))
		MR = H21@la.inv(E*np.eye(N)-H11)@H12@(Fp2+np.eye(N))

		return ML, MR

	def self_energy_RGF(self,E) :
		H11, H12, H21, H22 = self.hamiltonian()
		O = np.zeros(np.shape(H11))
		beta_dag = np.block([[O, H12],[O, O]])
		beta = np.block([[O, O],[H21, O]])
		alpha = np.block([[H11, H12],[H21, H22]])

		sz = np.shape(alpha)
		sz = sz[0]
		sub = np.shape(H11)
		sub = sub[0]

		eta = 1e-6

		g1_prev = np.zeros(np.shape(alpha))
		error = 100
		while error > 1e-8 :
			Gs = la.inv(complex(E,eta)*np.eye(sz)-alpha-(beta_dag@g1_prev@beta))
			error = np.sum(np.sum(np.abs(Gs-g1_prev)))/np.sum(np.sum(np.abs(Gs)+np.abs(g1_prev)))
			g1_prev = (Gs+g1_prev)*0.5

		M1 = beta_dag@g1_prev@beta

		g2_prev = np.zeros(np.shape(alpha))
		error = 100
		while error > 1e-8 :
			Gs = la.inv(complex(E,eta)*np.eye(sz)-alpha-(beta@g2_prev@beta_dag))
			error = np.sum(np.sum(np.abs(Gs-g2_prev)))/np.sum(np.sum(np.abs(Gs)+np.abs(g2_prev)))
			g2_prev = (Gs+g2_prev)*0.5

		M2 = beta@g2_prev@beta_dag

		return M1[:sub,:sub], M2[-sub:,-sub:]


	def self_energy_modes(self,E) :
		Up1, Up2, Un1, Un2, wp, wn = self.subcolumn_modes(E)
		H11, H12, H21, H22 = self.hamiltonian()
		N = len(wp)

		MLout_mode = []
		MRout_mode = []
		MLin_mode = []
		MRin_mode = []	

		Up1_inv = la.inv(Up1)
		Un1_inv = la.inv(Un1)
		Up2_inv = la.inv(Up2)
		Un2_inv = la.inv(Un2)

		for i in range(N) :
			MLout_mode.append( H12@la.inv(E*np.eye(N)-H22)@H21@(((1/wn[i])+1)*(Un1[:,i:i+1]*Un1_inv[i:i+1,:])) ) 
			MLin_mode.append( H12@la.inv(E*np.eye(N)-H22)@H21@((wp[i]+1)*(Up1[:,i:i+1]*Up1_inv[i:i+1,:])) )
			MRin_mode.append( H21@la.inv(E*np.eye(N)-H11)@H12@(((1/wn[i])+1)*(Un2[:,i:i+1]*Un2_inv[i:i+1,:])) ) 
			MRout_mode.append( H21@la.inv(E*np.eye(N)-H11)@H12@((wp[i]+1)*(Up2[:,i:i+1]*Up2_inv[i:i+1,:])) )

		MLout_kp = np.zeros((N,N))
		MRout_kp = np.zeros((N,N))
		MLin_kp = np.zeros((N,N))
		MRin_kp = np.zeros((N,N))
		MLout_kn = np.zeros((N,N))
		MRout_kn = np.zeros((N,N))
		MLin_kn = np.zeros((N,N))
		MRin_kn = np.zeros((N,N))		
		
		for i in range(N) :
			k = np.log(wp[i])/1j
			if abs(np.imag(k)) <= 1e-10 :
				if np.real(k) > 0 :
					MRout_kp = MRout_kp + MRout_mode[i]
					MLin_kp = MLin_kp + MLin_mode[i]
				else : 
					MRout_kn = MRout_kn + MRout_mode[i]
					MLin_kn = MLin_kn + MLin_mode[i]

		for i in range(N) :
			k = np.log(wn[i])/1j
			if abs(np.imag(k)) <= 1e-10 :
				if np.real(k) > 0 :
					MLout_kp = MLout_kp + MLout_mode[i]
					MRin_kp = MRin_kp + MRin_mode[i]
				else : 
					MLout_kn = MLout_kn + MLout_mode[i]
					MRin_kn = MRin_kn + MRin_mode[i]

		return MLout_kp, MLout_kn, MRout_kp, MRout_kn, MLin_kp, MLin_kn, MRin_kp, MRin_kn



class blg_channel :
	def __init__(self, L, W, site_eps_fn,U,t,tL) :
		self.U = U
		self.t = t
		self.tL = tL
		self.num_columns = int(L)
		self.W = W
		self.subcolumns = []
		for  i in range(self.num_columns) :
			self.subcolumns.append(blg_subcolumn1(i,W,site_eps_fn))
			self.subcolumns.append(blg_subcolumn2(i+0.5,W,site_eps_fn))

		self.sites = []	
		for i in range(2*self.num_columns) :
			self.sites = self.sites + self.subcolumns[i].sites

	def hamiltonian(self) :
		N = len(self.sites)
		H = np.zeros((N,N))
		
		for i in range(N) :
			for j in range(N) :
				if i == j :
					H[i,j] = self.sites[i].onsite(self.U) 		
				else :
					H[i,j] = self.sites[i].hopping(self.sites[j],[-self.t,-self.tL])

		return H

	def green_fun(self, lead, E, fL, fR) :
		if self.W != lead.W :
			raise Exception

		H = self.hamiltonian()
		ML, MR = lead.self_energy(E)
		H11, H12, H21, H22 = lead.hamiltonian()

		O = np.zeros(np.shape(H11))
		beta_dag = np.block([[O, H12],[O, O]])
		beta = np.block([[O, O],[H21, O]])
		alpha = np.block([[H11, H12],[H21, H22]])
		alphaL = np.block([[H11+ML, H12],[H21, H22]])
		alphaR = np.block([[H11, H12],[H21, H22+MR]])

		sub = np.shape(ML)
		ch = np.shape(H)
		ld = np.shape(alpha)
		sub = sub[0]
		ch = ch[0]
		ld = ld[0]
		
		

		HH = np.block([[alphaL,np.zeros((ld,ch)),np.zeros((ld,ld))],
			[np.zeros((ch,ld)),H,np.zeros((ch,ld))],
			[np.zeros((ld,ld)),np.zeros((ld,ch)),alphaR]])
		
		HH[:ld,ld:2*ld] = beta 
		HH[ld:2*ld,:ld] = beta_dag
		HH[-ld:,-2*ld:-ld] = beta_dag 
		HH[-2*ld:-ld,-ld:] = beta

		sz = np.shape(HH)
		sz = sz[0]
		GR = la.inv(E*np.eye(sz)-HH)
		GA = np.transpose(np.conjugate(GR))

		A = 1j*(GR-GA)

		TL = (0+0j)*np.ones(np.shape(HH))
		TR = (0+0j)*np.ones(np.shape(HH))
		TL[:sub,:sub] = 1j*(ML - np.transpose(np.conjugate(ML)))
		TR[-sub:,-sub:] = 1j*(MR - np.transpose(np.conjugate(MR)))

		Gn = GR@(TL*fL+TR*fR)@GA
		Gp = A-Gn

		return GR, GA, A, Gn, Gp, TL, TR 	



class blg_system :
	def __init__(self,L,W,site_eps_fn_channel,site_eps_fn_lead,U,t,tL) :
		self.channel = blg_channel(L,W,site_eps_fn_channel,U,t,tL)
		self.lead = blg_lead(W,site_eps_fn_lead,t,tL)

	def Rcurrent(self,E,fL,fR) :
		GR, GA, A, Gn, Gp, TL, TR = self.channel.green_fun(self.lead,E,fL,fR)
		MLout_kp, MLout_kn, MRout_kp, MRout_kn, MLin_kp, MLin_kn, MRin_kp, MRin_kn = self.lead.self_energy_modes(E)

		sub = np.shape(MRout_kp)
		sub = sub[0]
		TRout_kp = 1j*(MRout_kp - np.transpose(np.conjugate(MRout_kp)))
		TRin_kp = 1j*(MRin_kp - np.transpose(np.conjugate(MRin_kp)))
		TRout_kn = 1j*(MRout_kn - np.transpose(np.conjugate(MRout_kn)))
		TRin_kn = 1j*(MRin_kn - np.transpose(np.conjugate(MRin_kn)))

		IR_kp = np.trace((1-fR)*TRout_kp@Gn[-sub:,-sub:]-fR*TRin_kp@Gp[-sub:,-sub:])
		IR_kn = np.trace((1-fR)*TRout_kn@Gn[-sub:,-sub:]-fR*TRin_kn@Gp[-sub:,-sub:])

		return IR_kp, IR_kn, np.trace(TR@GR@TL@GA)		

	def Lcurrent(self,E,fL,fR) :
		GR, GA, A, Gn, Gp, TL, TR = self.channel.green_fun(self.lead,E,fL,fR)
		MLout_kp, MLout_kn, MRout_kp, MRout_kn, MLin_kp, MLin_kn, MRin_kp, MRin_kn = self.lead.self_energy_modes(E)

		sub = np.shape(MLout_kp)
		sub = sub[0]
		TLout_kp = 1j*(MLout_kp - np.transpose(np.conjugate(MLout_kp)))
		TLin_kp = 1j*(MLin_kp - np.transpose(np.conjugate(MLin_kp)))
		TLout_kn = 1j*(MLout_kn - np.transpose(np.conjugate(MLout_kn)))
		TLin_kn = 1j*(MLin_kn - np.transpose(np.conjugate(MLin_kn)))

		IL_kp = np.trace((1-fL)*TLout_kp@Gn[:sub,:sub]-fL*TLin_kp@Gp[:sub,:sub])
		IL_kn = np.trace((1-fL)*TLout_kn@Gn[:sub,:sub]-fL*TLin_kn@Gp[:sub,:sub])

		return IL_kp, IL_kn, np.trace(TL@GR@TR@GA)	


			

def site_eps_fn(name,x,y) :
	return 0

lead = blg_lead(20*sqrt(3),site_eps_fn,1,0.5)
lead.plot_bandstructure(np.linspace(-pi,pi,100))

syst = blg_system(1,20*sqrt(3),site_eps_fn,site_eps_fn,0,1,0.5)
print(syst.Rcurrent(0.1,1,0))