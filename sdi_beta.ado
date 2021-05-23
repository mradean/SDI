
cap mata: mata drop acarg()
cap mata: mata drop atarg()
cap mata: mata drop at_xac()
cap mata: mata drop oarg()
cap mata: mata drop ovarg()
cap mata: mata drop roarg()
cap mata: mata drop sdi1()
cap mata: mata drop sdi2()
cap mata: mata drop sdi3()
cap mata: mata drop sdi4()
cap mata: mata drop sdi5()
cap mata: mata drop xvarg()
cap mata: mata drop dsxarg()

mata

	void dsxarg(__s_id,__s_ind,__s_ord,__s_rb,__s_rV) { // deyx() | svar | xofi() argument

		rng = strtoreal(st_local("rng"))
		ac = strtoreal(st_local("nac"))
		at = strtoreal(st_local("nat"))
		rb = st_matrix("r(b)")		
		al = cols(rb)
		rV = st_matrix("r(V)")

		nx = strtoreal(st_local("nq"))

		if (at==1) {
			id = J(nx,1,1)
			if (ac==1) {
				if (nx==2 | rng) st_local("lab","_xofi")
				else st_local("lab","Nº: _xofi ")
				ind = strofreal(J(1,1,1::nx))
			}
			else {
				if (nx==2 | !rng) st_local("lab","_acr#_xofi")
				else st_local("lab","Nº: _acr#_xofi ") 					
				ind = strofreal(J(1,1,1::nx)):+"#1"
			}
		}
		else {	
			id = vec(J(nx,1,1..at))
			if (ac==1) {
				st_local("lab","Nº: _at#_xofi ")
				ind = strofreal(id):+"#":+strofreal(J(at,1,1::nx))
			}
			else {
				st_local("lab","Nº: _at#_acr#_xofi ")
				ind = strofreal(id):+"#":+strofreal(J(at,1,1::nx)):+"#1"
			}
		}

		if (__s_ord[1]!=.) {
			rb  = rb[__s_ord]
			rV  = rV[__s_ord,__s_ord]
		};

		__s_id = id
		__s_ind = ind
		__s_rb = rb
		__s_rV = rV
	}	

	void acarg(__s_acind,__s_acname,__s_acoz,__s_naso,string scalar acr) { // across argument
		__s_acoz = colmissing(st_matrix("r(at)")):==0 // across, one and zeros
		__s_naso = selectindex(__s_acoz) // non "asobserved" stats

		if (cols(__s_naso)>1) {
			if (st_global("r(xvars)")=="") {
				det = st_matrixcolstripe("r(at)")[__s_naso,2] // column names; remove eq. #
			}
			else det = tokens(st_global("r(xvars)"))' // deyx tokens
			det = subinstr(det,"b.",".",.) // remove base indicator
			ind = substr(det,1,strpos(det,"."):-1) // factor indicator
			det = substr(det,strpos(det,"."):+1,.) // variable names without the factor indicator
			if (rows(uniqrows(det))>1) {
				errprintf("only one variable allowed with {bf:across}\n")
				exit(198)
			}
			__s_acname = det[1]
			st_strscalar(acr,__s_acname+"=("+invtokens(ind')+")")
			__s_acind = strtoreal(ind')
		}
		else {
			if (st_local("aci")!="") {
				errprintf("factor '"+__s_acname+"' not found in list of covariates\n")
				exit(198)
			}
			__s_acname = st_matrixcolstripe("r(at)")[__s_naso,2] // column names; remove eq. #
		}
	}

	void atarg(__s_atlv,__s_atm,__s_atname,__s_ats,__s_atoz,__s_fvar,__s_rat,__s_selv) { // at() argument
	
		__s_atlv=__s_rat = __s_atm[__s_selv,] // at() values
		atr = rows(__s_atlv)
		st_local("nat",strofreal(atr))
		cna = st_matrixcolstripe("r(at)")[,2] // column names; remove eq. #
		cna = subinstr(cna,"b.",".",.) // remove base indicator
		stat = tokens(st_global("r(atstats1)"))	

		if (cols(selectindex(__s_fvar))) {
			fvls = ((stat:=="base") :+ (stat:=="value") :+ (stat:=="values")):*(__s_fvar:*__s_atoz) // factor variables at value(s)
			baln = (stat:=="asbalanced"):*(__s_fvar:*__s_atoz)

			if (any(fvls)) {
				ind = strtoreal(substr(cna,1,strpos(cna,"."):-1)) // factor indicator
				id = substr(cna,strpos(cna,"."):+1,.)
				info = panelsetup(id,1,2)
				info = info[selectindex(fvls[info[,1]]),]
				info = info,(info[,2]:-1),(info[,2]-info[,1])
				ri = rows(info)
				
				for (i=1; i<=ri; i++) {
					__s_atlv[|1,info[i,1]\atr,info[i,2]|] = __s_atlv[|1,info[i,1]\atr,info[i,2]|]:*ind[|info[i,1]\info[i,2]|]'
					__s_atlv[,info[i,2]] = rowsum(__s_atlv[|1,info[i,1]\atr,info[i,2]|])
					cna[info[i,2]] = substr(cna[info[i,2]],strpos(cna[info[i,2]],"."):+1,.)
					__s_atoz[|info[i,1]\info[i,3]|] = J(1,info[i,4],0)
				}
			}
			else {
				id = ""
			}
			
			if (any(baln)) {
				if (id[1]=="") id = substr(cna,strpos(cna,"."):+1,.) ;;
				info = panelsetup(id,1,2)
				info = info[selectindex(baln[info[,1]]),]
				info = info,(info[,2]:-1),(info[,2]-info[,1])
				ri = rows(info)
				
				for (i=1; i<=ri; i++) {
					cna[info[i,2]] = substr(cna[info[i,2]],strpos(cna[info[i,2]],"."):+1,.)
					__s_atoz[|info[i,1]\info[i,3]|] = J(1,info[i,4],0)
				}
				baln = baln:*__s_atoz // keep just the last entry per group
			};
		};

		vals = ((stat:=="value"):+(stat:=="values")):*__s_atoz // factor variables at value(s)
		gnlb = (substr(stat,1,1):=="("):*__s_atoz // stat labels for gen()
		if (any(gnlb)) stat[selectindex(gnlb)] = substr(stat[selectindex(gnlb)],2,strlen(stat[selectindex(gnlb)]):-2) ;;
		if (cols(selectindex(__s_fvar))) {
			stat[selectindex(baln)] = "(":+stat[selectindex(baln)]:+")"
			gnlb = gnlb:+baln
		};
		stlb = __s_atoz:-(vals:+gnlb) // stat labels except values, base, and gen)
		stat[selectindex(vals)] = J(1,sum(vals),"")
		stat[selectindex(stlb)] = " (":+stat[selectindex(stlb)]:+")"

		selc = selectindex(__s_atoz)
		if (any(gnlb)) {
			__s_atlv = strofreal(__s_atlv)
			__s_atlv[selectindex(gnlb)] = substr(__s_atlv[selectindex(gnlb)],3,.) // remover missing indicator, ie "."
			__s_atlv = colshape(__s_atlv[,selc],1)+J(atr,1,stat[selc]')
		}
		else __s_atlv = colshape(strofreal(__s_atlv[,selc]),1)+J(atr,1,stat[selc]')
		st_local("ato",strofreal(rows(__s_atlv))) // # of obs for ats
		__s_atname = J(atr,1,cna[selc])

		atc = cols(selc)
		if (atc==1 & atr==1) __s_ats = "at :"
		else if (atc>1 & atr==1) __s_ats = "at :"\J(atc-1,1,"")
		else if (atc==1 & atr>1) __s_ats = strofreal(1::atr):+"._at :"
		else {
			__s_ats = strofreal(1::atr):+"._at :",J(atr,atc-1,"")
			__s_ats = colshape(__s_ats,1)
		}
	}

	void at_xac(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac,__s_selv) { // when at() is used in with xval and/or across

		stat = tokens(st_global("r(atstats1)"))
		mxvs = max(selectindex(stat:=="values")) // maximum col # with valueS
		if (mxvs!=. & mxvs>max(__s_naso)) oarg(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac)
		else if (st_local("across")!="" & st_local("over")=="") {
			__s_ord = 0
			oarg(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac) 
		};
		
		if (st_local("`nolegend'")=="" & st_local("over")=="") {
			atr = rows(__s_atm)	
			if (cols(__s_naso)>1) {
				__s_selv = rowsum(__s_atm[,__s_naso]:*__s_acind),(1::atr) // use the numeric indicator for factor variables
			}
			else __s_selv = __s_atm[,__s_naso],(1::atr)
			fstr = __s_selv[,1]:==__s_selv[1,1] // selecting rows with the 1st value
			__s_selv = __s_selv[selectindex(fstr),2]
		}
		else if (st_local("`nolegend'")=="" & st_local("over")!="") {
			ac = strtoreal(st_local("nac"))
			__s_selv = range(1,(rows(__s_atm)/ac-1)*ac+1,ac) 
		};
	}

	void oarg(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac) { // order argument

		rng = strtoreal(st_local("rng"))
		atr = rows(__s_atm)
		if (cols(__s_naso)>1) group = rowsum(__s_atm[,__s_naso]:*__s_acind),(1::atr) // use the numeric indicator for factor variables
		else {

		group = __s_atm[,__s_naso],(1::atr) // group(across)
		}
		_sort(group,(1,2))
		sr = any(group[,2]:!=(1::atr)) // sorting required
		info = panelsetup(group,1)
		ngr = info[1,2] // # of elements in a group()
		ri = rows(info) // # of groups
			
		if (!rng & sr) {
			info = info,group[info[,1],2]
			_sort(info,3) // sort panels not by value, but by occurrence
		} ;
		if (st_local("across")!="" & st_local("over")=="") {			
			if (ri<2) {
				errprintf("need more than 2 for pairwise comparison\n")
				exit(198)
			};
			if (st_local("xno")!="1") {
				errprintf("only one {it:xofivar} value allowed with {bf:across}\n")
				exit(198)
			};
			
			st_local("nac",strofreal(ri)) // # of across()
			if (st_local("nolegend")=="") {
				__s_acat = strofreal(J(1,1,1::ri)):+"._acr :"
				__s_acname = J(ri,1,__s_acname)
				__s_rac = group[info[,1],1] // levels of across
				__s_aclv = strofreal(__s_rac)
			}
		};
		
		if (__s_ord | (!__s_ord & (rng & sr))) {
			group = group[,2]
			__s_ord = J(ngr,ri,.) // order vector
			for (i=1; i<=ri; i++) {
				__s_ord[|1,i\ngr,i|] = panelsubmatrix(group,i,info)
			}
			__s_ord = colshape(__s_ord,1)
		}
		else __s_ord = .
	}

	void ovarg(__s_acat,__s_aclv,__s_acname,__s_rac,string scalar es) { // over() argument
		st_view(lvs=.,.,__s_acname,es) // es = omit if !e(sample)
		__s_rac = uniqrows(lvs)
		__s_aclv = strofreal(__s_rac)
		ac = rows(__s_aclv)
		if (ac<2) {
			errprintf("need more than 2 for pairwise comparison\n")
			exit(198)
		};

		__s_acat = strofreal(1::ac):+J(ac,1,"._acr :")
		__s_acname = J(ac,1,__s_acname)
		st_local("nac",strofreal(ac))
	}

	void roarg(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac) { // order argument if rng==1 & acr>1 & at==1; to sort by acr

		rng = strtoreal(st_local("rng"))
		ac = rows(__s_atm)
		if (ac<2) {
			errprintf("need more than 2 for pairwise comparison\n")
			exit(198)
		};
		st_local("nac",strofreal(ac))

		if (st_local("nolegend")=="") {
			if (cols(__s_naso)==1) {
				if (!rng) __s_rac = __s_atm[,__s_naso]
				else {
					__s_aclv = sort((__s_atm[,__s_naso],(1::ac)),1)
					if (any(__s_aclv[,2]:!=(1::ac))) __s_ord = __s_aclv[,2] ;;
					__s_rac = __s_aclv[,1]
				}
			}
			else {
				if (!rng) __s_rac = rowsum(__s_atm[,__s_naso]:*__s_acind) // use the numeric indicator for factor variables
				else {
					__s_aclv = sort((rowsum(__s_atm[,__s_naso]:*__s_acind),(1::ac)),1)
					if (any(__s_aclv[,2]:!=(1::ac))) __s_ord = __s_aclv[,2] ;;
					__s_rac = __s_aclv[,1]
				}
			}
			__s_aclv = strofreal(__s_rac)
			__s_acat = strofreal(J(1,1,1::ac)):+"._acr :"
			__s_acname = J(ac,1,__s_acname)
		}
		else if (rng) {
			if (cols(__s_naso)>1) {
				__s_ord = sort((rowsum(__s_atm[,__s_naso]:*__s_acind),(1::ac)),1)
			}
			else {
				ord = sort((__s_atm[,__s_naso],(1::ac)),1)
			}
			if (all(__s_ord[,2]:==(1::ac))) __s_ord = .
			else __s_ord = __s_ord[,2]
		};
	}

	// sdi-list; (df=0 | dyx!=""), acr==0, range==0
	void sdi1(__s_df,__s_id,__s_ind,__s_lab,__s_matn,__s_rlab,__s_rb,__s_rV,__s_typ) {

		rlev = strtoreal(st_local("level"))/100
		if (__s_df==.) zs = invnormal(1-.5*(1-rlev))
		else zs = invttail(__s_df,(.5*(1-rlev)))
		prec = 1/10^(strtoreal(st_local("precision")))
		mc = strtoreal(st_local("mc"))
		rb = __s_rb'
		seb = sqrt(diagonal(__s_rV)) // b's standard errors 
	
		// creating the gradient for the first difference, gf
		rid = rows(__s_id)
		if (rid==2) {
			ri = 1
			pair = 1,2
			gf = (-1,1)
		}
		else {
			d = J(1,rid,1)
			g = diag(d)	
	
			info = panelsetup(__s_id,1,2)
			n  = info[1,2] // number of combination elements
			ri = rows(info)
			
			if (n==2) {
				k = 1
				pair = 1,2
			}
			else {
				n1 = n-1
				k = n*n1/2 // number of combinations
				pair = J(k,2,.)

				v = 1::n
				c1 = 1
				c2 = n1
				for (j=1; j<n; j++) {
					pair[|c1,1\c2,2|] = J(n-j,1,v[j]) , v[j+1..n]
					c1 = c2 +1
					c2 = c2 +n1 -j
				}
			}
			
			if (ri>1) {
				pair = J(ri,1,pair):+(vec(J(k,1,0..(ri-1)))*n)
			};

			gf = g[pair[,2],]:-g[pair[,1],] // gradient, first difference
		}
		vcf = gf*__s_rV*gf' // variance-covariance matrix, first diff		
		se_v = sqrt(diagonal(vcf)) // standard errors, first diff
		z_d = (zs:*se_v:+strtoreal(st_local("mvalue"))):/rowsum((seb[pair[,1]],seb[pair[,2]]))		
		if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
		else {
			st_local("cens","1")
			if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
			if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
		}

		rp = rows(pair)
		if (st_local("difference")=="") {			
			__s_rlab = 1::(2+2*rp)

			cpair = colshape(pair,1) // one column pairs
			qoi = rb[cpair]  // quantity of interest
			seq = seb[cpair] // first difference
			lab = __s_ind[cpair] // labels
			if (mc) lab = strofreal(vec(J(2,1,1..rp))):+": ":+lab ;;
			sdi = colshape(J(1,2,(1:-2*(1:-(normal(z_d))))*100),1)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(trunc(pdi):!=pdi)) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // trunc select index
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			typ = strofreal(pdi, st_local("form")):+J(2*rp,1," SDI")
			z_d = colshape((z_d,z_d),1)			
		}
		else {
			__s_rlab = 1::(2+3*rp)
		
			rb_v = rowsum(gf*rb) // first difference results
			qoi = colshape((rb[pair[,1]],rb[pair[,2]],rb_v),1)
			seq = colshape((seb[pair[,1]],seb[pair[,2]],se_v),1)
			if (st_local("deyx")=="") lab = (__s_ind[pair[,1]],__s_ind[pair[,2]]):+J(rp,1,(" (a)"," (b)"))
			else {
				__s_ind = strrtrim(__s_ind)
				lab = (__s_ind[pair[,1]],__s_ind[pair[,2]]):+J(rp,1,(" (a)  "," (b)  "))
			}
			if (!mc) lab = colshape((lab,J(rp,1,("(b)vs(a)"))),1)		 
			else {
				lab = colshape((lab,strofreal(1::rp):+".":+J(rp,1,("(b)vs(a)"))),1)		 
				lab = strofreal(vec(J(3,1,1..rp))):+": ":+lab
			}
			sdi = J(1,2,(1:-2*(1:-(normal(z_d))))*100)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(pdi:!=trunc(pdi))) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // select index for trunc
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			cil = J(rp,1,st_local("level")) // CI level	
			typ = colshape((strofreal(sdi,st_local("form")),cil),1):+J(rp,1,(" SDI"\" SDI"\"  CI"))
			sdi = colshape((sdi,strtoreal(cil)),1)
			z_d = colshape((z_d,z_d,J(rp,1,zs)),1)
		}

		pm = z_d:*seq // plus minus
		lo = qoi-pm
		hi = qoi+pm

		ztv = qoi:/seq // z- or t-values
		if (__s_df==.) pval = 2*(normal(-abs(ztv)))
		else pval = 2*ttail(__s_df,abs(ztv))
		__s_matn = J(2,8,.)\(qoi,seq,ztv,pval,z_d,lo,hi,sdi)
		__s_lab = st_local("lab")\" "\lab
		__s_typ = J(2,1,"")\typ		
	}
	
	// sdi-list; (pwc=0 | dyx!=""), acr==0, range==1
	void sdi2(__s_df,__s_id,__s_ind,__s_lab,__s_matn,__s_rlab,__s_rb,__s_rV,__s_typ) {		

		rlev = strtoreal(st_local("level"))/100
		if (__s_df==.) zs = invnormal(1-.5*(1-rlev))
		else zs = invttail(__s_df,(.5*(1-rlev)))
		prec = 1/10^(strtoreal(st_local("precision")))
		mc = strtoreal(st_local("mc"))
		rb = __s_rb'
		seb = sqrt(diagonal(__s_rV)) // b's standard errors 		
		
		rid = rows(__s_id) // # of id rows
		info = panelsetup(__s_id,1,2) // extreme values for panels
		n  = info[1,2] // number of combination elements
		ri = rows(info)
		
		inbw = info:+(1,-1) // the in-between values
		srt = (info[1,]')\(inbw[1,1]::inbw[1,2])
		gf = J(ri,rid,0) // gradient, first difference

		if (ri==1) {
			gf[info] = (-1,1)
		}
		else {
			srt = J(ri,1,srt):+(vec(J(n,1,0..(ri-1)))*n)
			for (i=1; i<=ri; i++) {
				gf[i,info[i,]] = (-1,1)
			}
		}

		vcf = gf*__s_rV*gf' // variance-covariance matrix, first diff
		se_v = sqrt(diagonal(vcf)) // standard errors, first diff
		
		if (st_local("difference")=="") {
			z_d = (zs:*se_v:+strtoreal(st_local("mvalue"))):/rowsum((seb[info[,1]],seb[info[,2]]))
			if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
			else {
				st_local("cens","1")
				if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
				if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
			}

			dv = colshape(info,1) // row #s of the data vector
			rv = colshape((info[,1],inbw[,1]),1) // row #s of the results vector
			seq=zr=lo=hi=sdir = J(rid,1,.)
			typ=lab = J(rid,1,"")

			qoi = rb[srt]
			seq[rv] = seb[dv]
			sdi = colshape(J(1,2,(1:-2*(1:-(normal(z_d))))*100),1)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(trunc(pdi):!=pdi)) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // trunc select index
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			typ[rv] = strofreal(pdi,st_local("form")):+J(ri,1,(" SDI"\" SDI"))
			z_d = colshape((z_d,z_d),1)
			lab = __s_ind[srt]
			if (mc) lab = strofreal(vec(J(n,1,1..ri))):+": ":+lab ;;
		}	
		else {
			rb_v = rowsum(gf*rb) // first difference results
			z_d = (zs:*se_v:+strtoreal(st_local("mvalue"))):/rowsum((seb[info[,1]],seb[info[,2]]))
			if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
			else {
				st_local("cens","1")
				if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
				if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
			}

			rid = rid+ri // add rows for contrast
			ri0 = 0::(ri-1)
			ri1 = 1::ri
			inst = inbw[,1]:+ri1 // # of inserted rows
			main = J(rid,1,1)
			main[inst] = J(ri,1,0)
			main = selectindex(main) // non inserted, main row #s
			dv = colshape((info,(info[ri,2]:+ri1)),1) // row #s for the data vector
			rv = colshape(((info[,1]:+ri0),(inbw[,1]:+ri0),inst),1) // row #s for the results vector
			qoi=seq=zr=lo=hi=sdir = J(rid,1,.)
			typ=lab = J(rid,1,"")

			qoi[main] = rb[srt]
			qoi[inst] = rb_v
			seq[rv] = (seb\se_v)[dv]
			sdi = J(1,2,(1:-2*(1:-(normal(z_d))))*100)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(pdi:!=trunc(pdi))) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // select index for trunc
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			cil = J(ri,1,st_local("level")) // CI level	
			typ[rv] = colshape((strofreal(sdi,st_local("form")),cil),1):+J(ri,1,(" SDI"\" SDI"\"  CI"))
			sdi = colshape((sdi,strtoreal(cil)),1)
			z_d = colshape((z_d,z_d,J(ri,1,zs)),1)
			if (st_local("deyx")=="") {
				lab[main] = __s_ind[srt]
				lab[rv] = lab[rv]:+J(ri,1,(" (a)"\" (b)"\"(b)vs(a)"))
			}
			else {
				lab[main] = strrtrim(__s_ind[srt])
				lab[rv] = lab[rv]:+J(ri,1,(" (a)  "\" (b)  "\"(b)vs(a)"))
			}
			if (mc) {
				lab[inst] = strofreal(ri1):+".(b)vs(a)"
				lab = strofreal(vec(J((n+1),1,ri1'))):+": ":+lab
			};
		}

		pm = z_d:*seq[rv] // plus minus
		lo[rv] = qoi[rv]-pm
		hi[rv] = qoi[rv]+pm

		ztv = qoi:/seq // z- or t-values
		if (__s_df==.) pval = 2*(normal(-abs(ztv)))
		else pval = 2*ttail(__s_df,abs(ztv))
		zr[rv] = z_d // z_d with range
		sdir[rv] = sdi
		__s_matn = J(2,8,.)\(qoi,seq,ztv,pval,zr,lo,hi,sdir)
		__s_lab = st_local("lab")\" "\lab
		__s_typ = J(2,1,"")\typ		
		__s_rlab = 1::(rid+2)
	}

	// NO FD
	void sdi3(__s_df,__s_lab,__s_matn,__s_ord,__s_rlab,__s_typ) {

		nx = strtoreal(st_local("nx"))
		rlev = strtoreal(st_local("level"))/100
		if (__s_df==.) zs = invnormal(1-.5*(1-rlev))
		else zs = invttail(__s_df,(.5*(1-rlev)))
		prec = 1/10^(strtoreal(st_local("precision")))
		mc = strtoreal(st_local("mc"))
		if (__s_ord[1]==.) {
			rb = st_matrix("r(b)")'
			rV = st_matrix("r(V)")
			rq = rows(rb)/2 // # of rows for quantity of interest
			at = rq/nx // # of at() values
		}
		else {
			rq = rows(__s_ord)
			at = rq/nx
			__s_ord = __s_ord\(__s_ord:+rq)
			rb = st_matrix("r(b)")'[__s_ord]		
			rV = st_matrix("r(V)")[__s_ord,__s_ord]
		}

		seb = sqrt(diagonal(rV)) // standard errors of the _b		
		gf  = diag(J(1,rq,-1)) , diag(J(1,rq,1)) // gradient, first difference
		rq1 = rq+1
		rq2 = rq*2
		qoi = rb[|1\rq|] , rb[|rq1\rq2|]
		seq = seb[|1\rq|] , seb[|rq1\rq2|]

		vcf = gf*rV*gf' // variance-covariance matrix, first diff
		sev = sqrt(diagonal(vcf)) // standard errors, first diff
		z_d = (zs:*sev:+strtoreal(st_local("mvalue"))):/rowsum(seq)
		if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
		else {
			st_local("cens","1")
			if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
			if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
		}

		if (st_local("difference")=="") {
			__s_rlab = 1::(2+rq2)

			qoi = colshape(qoi,1)
			seq = colshape(seq,1)
			sdi = colshape(J(1,2,(1:-2*(1:-(normal(z_d))))*100),1)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(trunc(pdi):!=pdi)) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // trunc select index
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			typ = strofreal(pdi,st_local("form")):+J(rq2,1," SDI")
			z_d = colshape((z_d,z_d),1)		
			num = strofreal(J(1,1,1::nx)) // numbering the pairs
			lab = num,("(":+num:+"+n)")
			if (at>1) {
				lab = J(at,1,lab)
				lab = strofreal(vec(J(2*nx,1,1..at))):+"#":+colshape(lab,1)
			}
			else lab = colshape(lab,1)
			if (mc) lab = strofreal(vec(J(2,1,1..at*nx))):+": ":+lab ;;
		}
		else {
			__s_rlab = 1::(2+3*rq)

			rb_v = rowsum(gf*rb) // second difference results
			qoi = colshape((qoi,rb_v),1)
			seq = colshape((seq,sev),1)
			sdi = J(1,2,(1:-2*(1:-(normal(z_d))))*100)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(pdi:!=trunc(pdi))) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // trunc select index
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			cil = J(rq,1,st_local("level")) // CI level	
			typ = colshape((strofreal(sdi,st_local("form")),cil),1):+J(rq,1,(" SDI"\" SDI"\"  CI"))
			sdi = colshape((sdi,strtoreal(cil)),1)
			z_d = colshape((z_d,z_d,J(rq,1,zs)),1)
			if (nx>1 | at>1) {
				num = strofreal(J(1,1,1::nx)) // numbering the pairs						
				lab = (num:+"     (a)"),("(":+num:+"+n) (b)"),(num:+".(b)vs(a)")
			}
			else lab = ("1     (a)"),("(1+n) (b)"),("(b)vs(a)")
			if (at>1) {
				lab = J(at,1,lab)
				lab = strofreal(vec(J(3*nx,1,1..at))):+"#":+colshape(lab,1)
			}
			else lab = colshape(lab,1)
			if (mc) lab = strofreal(vec(J(3,1,1..rq))):+": ":+lab ;;
		}
		
		pm  = z_d:*seq // plus minus
		lo  = qoi-pm
		hi  = qoi+pm

		ztv = qoi:/seq // z- or t-values
		if (__s_df==.) pval = 2*(normal(-abs(ztv)))
		else pval = 2*ttail(__s_df,abs(ztv))
		__s_matn = J(2,8,.)\(qoi,seq,ztv,pval,z_d,lo,hi,sdi)
		__s_lab = st_local("lab")\" "\lab		
		__s_typ = ""\""\typ
	}

	// FD (all)
	void sdi4(__s_df,__s_lab,__s_matn,__s_ord,__s_rlab,__s_typ) {

		nx = strtoreal(st_local("nx"))
		rlev = strtoreal(st_local("level"))/100
		if (__s_df==.) zs = invnormal(1-.5*(1-rlev))
		else zs = invttail(__s_df,(.5*(1-rlev)))
		prec = 1/10^(strtoreal(st_local("precision")))
		mc = strtoreal(st_local("mc"))
		if (__s_ord[1]==.) {
			rb = st_matrix("r(b)")'		
			rV = st_matrix("r(V)")
			rq = rows(rb)/2 // # of rows for quantity of interest
			at = rq/nx 
		}
		else {
			rq = rows(__s_ord)
			at = rq/nx
			__s_ord = __s_ord\(__s_ord:+rq)
			rb = st_matrix("r(b)")'[__s_ord]	
			rV = st_matrix("r(V)")[__s_ord,__s_ord]
		}
		
		gf = diag(J(1,rq,-1)) , diag(J(1,rq,1)) // gradient, first difference
		rb_v = rowsum(gf*rb) // first difference results
		vcf = gf*rV*gf' // variance-covariance matrix, first diff
		se_v = sqrt(diagonal(vcf)) // standard errors, first diff
		v = 1::nx
		if (st_local("nac")=="1") {
			if (at==1) ind = "[(":+strofreal(v):+"+n) vs ":+strofreal(v):+"]"
			else ind = strofreal(vec(J(nx,1,1..at))):+J(at,1,"#[(":+strofreal(v):+"+n) vs ":+strofreal(v):+"]")
		}
		else {
			if (at==1) ind = strofreal(v):+J(rq,1,"#[(1+n) vs 1]")
			else ind = strofreal(vec(J(nx,1,1..at))):+"#":+J(at,1,strofreal(v)):+J(rq,1,"#[(1+n) vs 1]")
		}
	
		if (nx==2 & at==1) {
			pair = 1,2
			gs = (-1,1) // gradient, second-first difference (treating fd as original results)
		}
		else {
			g = diag(J(1,rq,1))	
			n = nx // number of combination elements
			
			if (nx==2) {
				k = 1
				z = 1\2
				pair = 1,2
			}
			else {
				nx1 = nx-1
				k = nx*nx1/2 // number of combinations
				pair = J(k,2,.)

				c1 = 1
				c2 = nx1
				for (j=1; j<n; j++) {
					pair[|c1,1\c2,2|] = J(nx-j,1,v[j]) , v[j+1..nx]
					c1 = c2 +1
					c2 = c2 +nx1 -j
				}
			}
			
			if (at>1) {
				pair = J(at,1,pair):+(vec(J(k,1,0..(at-1)))*nx)
			};

			gs = g[pair[,2],]:-g[pair[,1],]
		}

		vcs = gs*vcf*gs' // variance-covariance matrix, second diff
		ses = sqrt(diagonal(vcs)) // standard errors, second diff
		z_d = (zs:*ses:+strtoreal(st_local("mvalue"))):/rowsum((se_v[pair[,1]],se_v[pair[,2]]))
		if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
		else {
			st_local("cens","1")
			if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
			if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
		}

		rp = rows(pair)
		if (st_local("difference")=="") {			
			__s_rlab = 1::(2+2*rp)
			
			cpair = colshape(pair,1) // one column pairs
			qoi = rb_v[cpair]  // quantity of interest
			seq = se_v[cpair] // first difference
			lab = ind[cpair] // labels
			if (mc) lab = strofreal(vec(J(2,1,1..rp))):+": ":+lab ;;
			sdi = colshape(J(1,2,(1:-2*(1:-(normal(z_d))))*100),1)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(trunc(pdi):!=pdi)) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // trunc select index
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			typ = strofreal(pdi,st_local("form")):+J(2*rp,1," SDI")
			z_d = colshape((z_d,z_d),1)			
		}
		else {
			__s_rlab = 1::(2+3*rp)

			rbs = rowsum(gs*rb_v) // second difference results
			qoi = colshape((rb_v[pair[,1]],rb_v[pair[,2]],rbs),1)
			seq = colshape((se_v[pair[,1]],se_v[pair[,2]],ses),1)
			lab = (ind[pair[,1]],ind[pair[,2]]):+J(rp,1,(" (a)"," (b)"))
			if (!mc) lab = colshape((lab,J(rp,1,("(b)vs(a)"))),1)		 		 
			else {
				lab = colshape((lab,strofreal(1::rp):+".":+J(rp,1,("(b)vs(a)"))),1)		 
				lab = strofreal(vec(J(3,1,1..rp))):+": ":+lab ;;
			}
			sdi = J(1,2,(1:-2*(1:-(normal(z_d))))*100)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(pdi:!=trunc(pdi))) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // select index for trunc
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			cil = J(rp,1,st_local("level")) // CI level	
			typ = colshape((strofreal(sdi,st_local("form")),cil),1):+J(rp,1,(" SDI"\" SDI"\"  CI"))
			sdi = colshape((sdi,strtoreal(cil)),1)
			z_d = colshape((z_d,z_d,J(rp,1,zs)),1)
		}

		pm = z_d:*seq // plus minus
		lo = qoi-pm
		hi = qoi+pm

		ztv = qoi:/seq // z- or t-values
		if (__s_df==.) pval = 2*(normal(-abs(ztv)))
		else pval = 2*ttail(__s_df,abs(ztv))
		__s_matn = J(2,8,.)\(qoi,seq,ztv,pval,z_d,lo,hi,sdi)
		__s_lab = st_local("lab")\" "\lab
		__s_typ = J(2,1,"")\typ		
	}

	// FD (range)
	void sdi5(__s_df,__s_lab,__s_matn,__s_ord,__s_rlab,__s_typ) {
		
		nx = strtoreal(st_local("nx"))
		rlev = strtoreal(st_local("level"))/100
		if (__s_df==.) zs = invnormal(1-.5*(1-rlev))
		else zs = invttail(__s_df,(.5*(1-rlev)))
		prec = 1/10^(strtoreal(st_local("precision")))
		mc = strtoreal(st_local("mc"))
		if (__s_ord[1]==.) {
			rb = st_matrix("r(b)")'		
			rV = st_matrix("r(V)")
			rq = rows(rb)/2 // # of rows for quantity of interest
			at = rq/nx
		}
		else {
			rq = rows(__s_ord)
			at = rq/nx // # of x or ac values
			__s_ord = __s_ord\(__s_ord:+rq)
			rb = st_matrix("r(b)")'[__s_ord]	
			rV = st_matrix("r(V)")[__s_ord,__s_ord]
		}

		gf = diag(J(1,rq,-1)) , diag(J(1,rq,1)) // gradient, first difference
		rb_v = rowsum(gf*rb) // first difference results		
		vcf = gf*rV*gf' // variance-covariance matrix, first diff
		se_v = sqrt(diagonal(vcf)) // standard errors, first diff
		v = 1::nx
		if (st_local("nac")=="1") {
			if (at==1) ind = "[(":+strofreal(v):+"+n) vs ":+strofreal(v):+"]"
			else ind = strofreal(vec(J(nx,1,1..at))):+J(at,1,"#[(":+strofreal(v):+"+n) vs ":+strofreal(v):+"]")
		}
		else {
			if (at==1) ind = strofreal(v):+J(rq,1,"#[(1+n) vs 1]")
			else ind = strofreal(vec(J(nx,1,1..at))):+"#":+J(at,1,strofreal(v)):+J(rq,1,"#[(1+n) vs 1]")
		}
		
		info = 1,nx
		inbw = info:+(1,-1) // the in-between values
		srt = info'\(inbw[1]::inbw[2])
		gs = J(at,rq,0)

		if (at==1) {
			gs[info] = (-1,1) // gradient, second-first difference (treating fd as origianal results)
		}
		else {
			srt = J(at,1,srt):+vec(J(nx,1,0..(at-1)))*nx
			v = (0::(at-1))*nx
			info = J(at,1,info):+v
			inbw = J(at,1,inbw):+v
			for (i=1; i<=at; i++) {
				gs[i,info[i,]] = (-1,1)
			}
		}

		vcs = gs*vcf*gs' // variance-covariance matrix, first-second diff
		ses = sqrt(diagonal(vcs)) // standard errors, first-second diff
		
		if (st_local("difference")=="") {
			z_d = (zs:*ses:+strtoreal(st_local("mvalue"))):/rowsum((se_v[info[,1]],se_v[info[,2]]))
			if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
			else {
				st_local("cens","1")
				if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
				if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
			}

			dv = colshape(info,1) // row #s of the data vector
			rv = colshape((info[,1],inbw[,1]),1) // row #s of the results vector
			seq=zr=lo=hi=sdir = J(rq,1,.)
			typ=lab = J(rq,1,"")

			qoi = rb_v[srt]
			seq[rv] = se_v[dv]
			sdi = colshape(J(1,2,(1:-2*(1:-(normal(z_d))))*100),1)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(trunc(pdi):!=pdi)) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // trunc select index
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			typ[rv] = strofreal(pdi,st_local("form")):+J(at,1,(" SDI"\" SDI"))
			z_d = colshape((z_d,z_d),1)
			lab = ind[srt]
			if (mc) lab = strofreal(vec(J(nx,1,1..at))):+": ":+lab ;;
		}	
		else {
			rbs = rowsum(gs*rb_v) // second difference results
			z_d = (zs:*ses:+strtoreal(st_local("mvalue"))):/rowsum((se_v[info[,1]],se_v[info[,2]]))
			if (min(z_d)>=0 & max(z_d)<=5.85) st_local("cens","0") // censored z-score
			else {
				st_local("cens","1")
				if (min(z_d)<0) z_d = (z_d:>0) :* z_d ;;
				if (max(z_d)>5.85) z_d = ((z_d:>5.85)*5.85) :+ ((z_d:<5.85):*z_d) ;;
			}

			rq = rq+at // add rows for df
			at0 = 0::(at-1)
			at1 = 1::at
			inst = inbw[,1]:+at1 // # of inserted rows
			main = J(rq,1,1)
			main[inst] = J(at,1,0)
			main = selectindex(main) // non inserted, main row #s
			dv = colshape((info,(info[at,2]:+at1)),1) // row #s for the data vector
			rv = colshape(((info[,1]:+at0),(inbw[,1]:+at0),inst),1) // row #s for the results vector
			qoi=seq=zr=lo=hi=sdir = J(rq,1,.)
			typ=lab = J(rq,1,"")

			qoi[main] = rb_v[srt]
			qoi[inst] = rbs
			seq[rv] = (se_v\ses)[dv]
			sdi = J(1,2,(1:-2*(1:-(normal(z_d))))*100)
			pdi = sdi/prec // precision-corrected difference interval
			if (all(pdi:!=trunc(pdi))) pdi = prec*floor(pdi):+prec
			else {
				trcs = (trunc(pdi):!=pdi) // select index for trunc
				pdi[selectindex(trcs)] = prec*floor(pdi[selectindex(trcs)]):+prec
				pdi[selectindex(1:-trcs)] = prec*floor(pdi[selectindex(1:-trcs)])				
			}
			cil = J(at,1,st_local("level")) // CI level		
			typ[rv] = colshape((strofreal(sdi,st_local("form")),cil),1):+J(at,1,(" SDI"\" SDI"\"  CI"))
			sdi = colshape((sdi,strtoreal(cil)),1)
			z_d = colshape((z_d,z_d,J(at,1,zs)),1)
			lab[main] = ind[srt]
			lab[rv] = lab[rv]:+J(at,1,(" (a)"\" (b)"\"(b)vs(a)"))
			if (mc) {
				lab[inst] = strofreal(at1):+".(b)vs(a)"
				lab = strofreal(vec(J((nx+1),1,at1'))):+": ":+lab
			};
		}
		
		pm = z_d:*seq[rv] // plus minus
		lo[rv] = qoi[rv]-pm
		hi[rv] = qoi[rv]+pm

		ztv = qoi:/seq // z- or t-values
		if (__s_df==.) pval = 2*(normal(-abs(ztv)))
		else pval = 2*ttail(__s_df,abs(ztv))
		zr[rv] = z_d // z_d with range
		sdir[rv] = sdi
		__s_matn = J(2,8,.)\(qoi,seq,ztv,pval,zr,lo,hi,sdir)
		__s_lab = st_local("lab")\" "\lab
		__s_typ = J(2,1,"")\typ
		__s_rlab = 1::(rq+2)
	}

	void xvarg(__s_acoz,__s_atoz,__s_naso,__s_ovoz,__s_rx,__s_xat,__s_xlev,__s_xname,__s_xvoz,string scalar es,string scalar xvl1,string scalar xvl2) {
	
		nl = (st_local("nolegend")=="")
		nu = (st_local("nunit")!="")
		if (st_global("r(atstats1)")=="") { // currently (asobserved); run margins, at((mean) xofivar) to get __s_xname
			st_local("atx",ustrregexra(st_local("xofinterest"),"\(.+?\)",""))
			stata("margins in "+st_local("fsamp")+", at("+st_local("atx")+") nose noatlegend")
			cna = st_matrixcolstripe("r(at)")[.,2] // column names; remove eq. #
			__s_xvoz = colmissing(st_matrix("r(at)")):==0 // xofivar, one and zeros
			xvnao = selectindex(__s_xvoz) // non "asobserved" cols
			if (cols(xvnao)>1) {
				errprintf("only one non-factor variable allowed with {bf:xofinterest}\n")
				exit(198)
			};
			__s_xname = cna[xvnao] // name of xofivar
			if (st_local("nunit")!="") {
				if (st_local("sd")=="3") unit = strtoreal(st_local("nunit"))
				else {
					stata("sum "+__s_xname+" if "+es)
					if (st_local("sd")=="1") unit = st_numscalar("r(sd)")
					else unit = -st_numscalar("r(sd)")				
				}
				st_local("unit",strofreal(unit))
			};
			xno = 1 // number of xofivar values

			st_strscalar(xvl1,"(asobserved)"+__s_xname)
			if (nu) st_strscalar(xvl2,__s_xname+"=gen("+__s_xname+"+"+strofreal(unit)+")") ;;
			if (nl) {
				__s_rx = .o
				__s_xlev = "(asobserved)" 
			};
		}
		else {	
			cna = st_matrixcolstripe("r(at)")[.,2] // column names; remove eq. #
			atm = st_matrix("r(at)")			
			stat = tokens(st_global("r(atstats1)"))
			__s_xvoz = colmissing(atm[1,]):==0 // xofivar, one and zeros
			atm:!=max(atm)
			xvnao = selectindex(__s_xvoz)	// non "asobserved" stats
			if (cols(xvnao)>1) {
				errprintf("only one non-factor variable allowed with {bf:xofinterest}\n")
				exit(198)
			}
			else if (!cols(xvnao)) {
				__s_xvoz = atm:!=max(atm) // for the gen() case with a '.g' value
				xvnao = selectindex(__s_xvoz)
			};
			__s_xname = cna[xvnao] // name of xofivar

			if (st_local("nunit")!="") {
				if (st_local("sd")=="3") unit = strtoreal(st_local("nunit"))
				else stata("sum "+__s_xname+" if "+es)
				if (st_local("sd")=="1") unit = st_numscalar("r(sd)")
				else if (st_local("sd")=="2") unit = -st_numscalar("r(sd)")	;;
				st_local("unit",strofreal(unit))
			};
			xstat = stat[xvnao]
			xno = rows(atm) // number of xofivar values
			
			if (xstat=="value" | xstat=="values") {
				__s_xlev = atm[,xvnao]
				if (st_local("rng")=="1") _sort(__s_xlev,1) ;;
				st_strscalar(xvl1,__s_xname+"=("+invtokens(strofreal(__s_xlev'))+")")
				if (nu) {
					__s_xlev2 = __s_xlev:+unit
					st_strscalar(xvl2,__s_xname+"=("+invtokens(strofreal(__s_xlev2'))+")")
				};
				if (nl) {
					__s_rx = __s_xlev
					__s_xlev = strofreal(__s_xlev)
				};
			}
			else if (substr(xstat,1,1)!="(") {
				if (xstat=="median") xstat = "p50" ;;
				if (substr(xstat,1,1)=="m") {
					stata("sum "+__s_xname+" if "+es+", meanonly")
					__s_xlev = st_numscalar("r("+xstat+")")
				}
				else if (substr(xstat,1,1)=="p") {
					stata("_pctile "+__s_xname+" if "+es+", p("+xstat+")")
					__s_xlev = st_numscalar("r(r1)")
				}
				else {
					__s_xlev = 0
				}

				st_strscalar(xvl1,__s_xname+"=("+strofreal(__s_xlev)+")")
				if (nu) {
					__s_xlev2 = __s_xlev:+unit				
					st_strscalar(xvl2,__s_xname+"=("+strofreal(__s_xlev2)+")")
				};
				if (nl) {
					__s_rx = __s_xlev					
					__s_xlev = strofreal(__s_xlev)+" ("+stat[xvnao]+")"	
				};
			}
			else {
				__s_xlev = xstat
				st_strscalar(xvl1,__s_xname+"=gen"+__s_xlev)
				if (nu) st_strscalar(xvl2,__s_xname+"=gen("+__s_xlev+"+"+strofreal(unit)+")") ;;
				if (nl) __s_rx = .g ;;
			}
			
			if (nl) __s_xname = J(xno,1,__s_xname) ;;
		}
	
		if (nu) {
			if (xno==1 & __s_acoz[1]==. & __s_ovoz[1]==. & st_local("fd")=="1") {
				errprintf("need more than 2 for pairwise comparison\n")
				exit(198)
			};
		}
		else {
			if (xno==1 & __s_acoz[1]==. & __s_ovoz[1]==.) {
				errprintf("need more than 2 for pairwise comparison\n")
				exit(198)
			};
		}

		st_local("xno",strofreal(xno))
		if (nl) __s_xat = strofreal(J(1,1,1::xno)):+"._xofi :" ;;
		
		if (any((__s_acoz[1],__s_atoz[1],__s_ovoz[1]):!=.)) {
			if (__s_acoz[1]==.) {
				__s_acoz = J(1,cols(__s_xvoz),0)
				__s_naso = xvnao
			};
			if (__s_atoz[1]==.) __s_atoz = J(1,cols(__s_xvoz),0) ;;
			if (__s_ovoz[1]==.) __s_ovoz = J(1,cols(__s_xvoz),0) ;;
			if (max(__s_xvoz+__s_acoz+__s_atoz+__s_ovoz)>1) {
				errprintf("variables in {bf:xofinterest}, {bf:across}, and {bf:at} must be distinct\n")
				exit(198)
			};
		}
		else if (__s_acoz[1]==. & __s_ovoz[1]==.) __s_naso = xvnao ;;
	}
end


cap program drop sdi
program sdi
//version 15

	syntax [if] [in] [fw aw pw iw] , XOFInterest(str) [ACRoss(str) ASBALanced atmeans at(str) ///
			CONTinuous DIFFerence EXPression(str) fd firstdiff /// 
			Level(cilevel) many Mvalue(real 0) Nolegend Nonotes NOWEIGHTs nunit(str) PRedict(str) /// 
			PRECision(integer 1) range]

	cap noi {
		qui {
			local cmdl `"sdi `*'"'

			cap mata: st_local("sname",strofreal(max(substr(direxternal("*"),1,4):=="__s_"))) 
			if (!_rc & "`sname'"!="0") {
				di as err "There are Mata matrixes that start with '__s_' that you must either drop or rename before proceeding." ///
				_n "Use either {it:mata drop} or {it:mata rename}; see {ul:[M-3] Commands for controlling Mata}." ///
				_n "(To drop all, you can write {cmd:mata: mata drop __s_*} in Stata's command window.)" //?? SMCL
				exit(198)
			}

			// check that the previous command was not "margins, post"	
			if "`e(cmd)'"=="margins" {
				di as err "{bf:sdi} cannot work with {bf:margins}' posted results"
				exit(198)
			}
			
			// 'options' argument
			local options noatlegend predict(`predict') `asbalanced' `atmeans' `continuous' `noweights'
			
			local N _N
			estimates esample
			if r(who)=="zero'd" {
				di as err "{bf:e(sample)} does not identify the estimation sample"
				exit(198)
			}
			tempvar obsno esamp
			gen long `obsno' = _n
			gen long `esamp' = 0
			if ("`if'"!="") replace `esamp' `in' `if' & e(sample) = 1
			else replace `esamp' `in' if e(sample) = 1
			cap sum `obsno' if `esamp', meanonly
			if ("`r(min)'"=="") {
				di as err "no observations"
				exit(198)
			}
			else local fsamp `r(min)' // first observation in sample
			margins in `fsamp', predict(`predict') nose
				if (colsof(r(b))>1) {
					di as err "models with multiple outcomes and multi-equation models are not allowed"
					exit(198)
				}
			
			if "`range'"=="" local rng 0
			else local rng 1

			if ("`firstdiff'"!="" | "`fd'"!="") local fd 1
			else local fd 0

			tempname acro
			scalar `acro' = ""
			mata: __s_acat=__s_acind=__s_aclv=__s_acname=__s_acoz=__s_rac=__s_naso=__s_ovoz = . // ac is acno: # of across levels; sdiX() are written using ac
			if "`across'"!="" {
				if (strpos("`across'","#")!=0 ) {
					di as err "invalid {bf:across} option; levels of interactions not allowed"
					exit(198)
				}

				if (strpos("`across'","asobserved")) {	
					local acr = subinstr(stritrim("`across'"), "= ", "=",.)
					local acr = subinstr("`acr'","( ","(",.)
					local acr = subinstr("`acr'"," )",")",.)
					
					local paso = strpos("`acr'","(asobserved)") // position of (asobs)
					local pgen = strrpos("`acr'","=gen")
					if !(`pgen') {
						di as err "in {bf:across}, suboption {bf:(asobserved)} requires a matching {bf:generate()} suboption"
						exit(198)
					}
					else {
						if (`paso'<`pgen') {
							local acr = subinstr("`acr'", " =gen", "=gen",.)
							local pgen3 = strrpos(substr("`acr'",1,`pgen')," ")	// space position that separates gen()-suboption	
							if (!`pgen3') local pgen3 = strrpos(substr("`acr'",1,`pgen'),")")	
							if (`pgen3'==12) local pgen3 13
							margins in `fsamp', at(`=substr("`acr'",`paso'+12,`pgen3'-`paso'-12)') at(`=substr("`acr'",`pgen3',.)') nose noatlegend // name of asobs var
							if (`pgen3'==13) { // it means there is no varname after (asobserved)
								di as err "need more than 2 for pairwise comparison"
								exit(198)
							}
							// mata: acarg() by hand
							mata: __s_acoz = (st_matrix("r(at)"):!=.)
							mata: __s_acoz = colmax(__s_acoz) // across, one and zeros
							mata: __s_acname = J(2,1,st_matrixcolstripe("r(at)")[selectindex(__s_acoz),2]) // column names; remove eq. #
							mata: st_local("nacnm",strofreal(rows(__s_acname))) // number of (var) names
							if (`nacnm'>2) {
								di as err "only one variable allowed with {bf:across}"
								exit(198)
							}
							mata: __s_acat = "1._acr : "\"2._acr : "
							local gval = substr("`acr'",`pgen',.) // value of gen()
							local pgen1 `=strpos("`gval'","(")+1'
							local pgen2 `=strrpos("`gval'",")")-`pgen1''
							local gval = substr("`gval'",`pgen1',`pgen2')
							mata: __s_aclv = "(asobserved)"\st_local("gval")
							mata: __s_rac = .o\.g
							
							local acasob1 "(asobserved) _all"						
							local acasob2 "`across'"						
						}
						else {
							cap margins in `fsamp', at(`=substr("`acr'",1,`paso'-1)') at(`=substr("`acr'",`paso'+12,.)') nose noatlegend
							if (_rc==198 & strlen(substr("`acr'",`paso',.))==12) { // it means there is no varname after (asobserved)
								di as err "need more than 2 for pairwise comparison"
								exit(198)
							}
							else if (_rc) margins in `fsamp', at(`=substr("`acr'",1,`paso'-1)') at(`=substr("`acr'",`paso'+12,.)') nose noatlegend
							// mata: acarg() by hand
							mata: __s_acoz = (st_matrix("r(at)"):!=.)
							mata: __s_acoz = colmax(__s_acoz) // across, one and zeros
							mata: __s_acname = J(2,1,st_matrixcolstripe("r(at)")[selectindex(__s_acoz),2]) // column names; remove eq. #
							mata: st_local("nacnm",strofreal(rows(__s_acname))) 
							if (`nacnm'>2) {
								di as err "only one variable allowed with {bf:across}"
								exit(198)
							}
							mata: __s_acat = "1._acr : "\"2._acr : "
							local gval = substr("`acr'",1,`paso'-1) // value of gen()
							local pgen1 `=strpos("`gval'","(")+1'
							local pgen2 `=strrpos("`gval'",")")-`pgen1''
							local gval = substr("`gval'",`pgen1',`pgen2')
							mata: __s_aclv = st_local("gval")\"(asobserved)"
							mata: __s_rac = .g\.o
							
							local acasob1 = substr("`acr'",1,`paso'-1)						
							local acasob2 "(asobserved) _all"						
						}
					}
					local nac 2
				}
				else {
					if (strpos(" `across'"," i.")) {
						local aci = subinstr(" `across'"," i."," ",.)					
						cap margins in `fsamp', at((base) `aci') dydx(`aci') nose noatlegend				
						if (_rc==198) margins `aci' in `fsamp', nose				
						else if (_rc) margins in `fsamp', at((base) `aci') dydx(`aci') nose noatlegend				
						
						mata: acarg(__s_acind,__s_acname,__s_acoz,__s_naso,"`acro'")
					}
					else if (strpos(" `across'"," o.")) {
						if ("`over'"!="") {
							di as err "option {bf:over()} not allowed"
							exit(198)
						}
						local over = subinstr(" `across'"," o."," ",.)
					}
					else {
						margins in `fsamp', at(`across') nose noatlegend
						mata: acarg(__s_acind,__s_acname,__s_acoz,__s_naso,"`acro'")
						scalar `acro' = "`across'"
					}
					local nac 2 // 2 instead of 1 indicates across is set; exact value set later		
				}				
			}
			else local nac 1

			if ("`over'"!="") {
				if (strpos("`over'","(")) local over `=subinstr("`over'","(","",.)'
				cap margins in `fsamp', over(`over') at(`over') nose noatlegend
				if (!_rc) {
					if (`: word count `r(over)''>1) {
						di as err "only one variable allowed with {bf:across}"
						exit(198)
					}
					else {
						mata: __s_acname=rover = st_global("r(over)")
						mata: __s_ovoz = colmissing(st_matrix("r(at)")):==0 // across, one and zeros
					}
				}
				else {
					margins in `fsamp', over(`over') nose noatlegend
					if (`: word count `r(over)''>1) {
						di as err "only one variable allowed with {bf:across}"
						exit(198)
					}
					else mata: __s_acname=rover = st_global("r(over)")
				}
			}

			if "`at'"!="" {
				if strpos("`at'","(")!=0 {
					local ata `=subinstr(stritrim("`at'"), "= ", "=",.)'
					local ata `=subinstr("`ata'", "=gen", "(=gen",.)'
					local ata `=subinstr("`ata'", "=(", "(=",.)'
					local ata `=ustrregexra("`ata'","\(.+?\)","")'
				}
				else local ata `at' // at as at(aspect)
				if ("`ata'"=="") {
					di as err "`at' invalid statistic"
					exit(198)
				}

				margins in `fsamp', at(`ata') at(_factor) nose noatlegend
				mata: __s_atoz = colmissing(st_matrix("r(at)")[1,]):==0 // at, ones and zeros
				mata: __s_fvar = colmissing(st_matrix("r(at)")[2,]):==0 // all factor variables: ones and zeros

				if ("`nolegend'"=="") mata: atlv=atname=ats = .
			}
			else {
				mata: __s_atoz = .
				local nat 1
			}
			
			mata: __s_xvoz=__s_svoz=__s_deoz = .
			if ("`xofinterest'"!="") {
				if (strpos("`xofinterest'","#")!=0 ) {
					di as err "invalid {bf:xofinterest} option; levels of interactions not allowed"
					exit(198)
				}

				if ("`nunit'"=="") {
					margins in `fsamp', at(`xofinterest') nose noatlegend
					tempname xval1 xval2
					mata: __s_rx=__s_xat=__s_xlev=__s_xname = .		
					mata: xvarg(__s_acoz,__s_atoz,__s_naso,__s_ovoz,__s_rx,__s_xat,__s_xlev,__s_xname,__s_xvoz,"`esamp'","`xval1'","`xval2'")
				}
				else if (!`fd' & `rng' & `nac'==1) {
					di as err "option {bf:range} not allowed"
					exit(198)
				}
				else mata: __s_xvoz = 0 // indicates xofivar() is not missing
			}
			
			if (`precision'<0 | `precision'>6) {
				di as err "{bf:precision} invalid -- invalid number, outside of allowed range"
				exit(198)
			}
			else local form = "%`=`precision'+4'.`precision'f"

			mata: __s_lab=__s_matn=__s_ord=__s_rlab=__s_typ = .

			// Results
			if ("`nunit'"=="") { // varlist!="" |  dydx!="" |  (xofi!="" & nunit=="")

				margins `varlist' `if' `in' `fw' `aw' `pw' `iw', `deyx' at(`=`acro'' `at' `=`xval1'') level(`level') over(`over') `options'		
				local df = (el(r(table),7,1)!=.)
				mata: __s_df = st_matrix("r(table)")[7,1]
				mata: __s_rN = st_numscalar("r(N)")
				if ("`nolegend'"=="") {
					local len `=strlen("r(N)")'
					if (`len'>3) local rno `=string(r(N),"%-"+"`=`len'+int(`len'/3)+2'"+".0gc")'
					else local rno `=string(r(N))'
					if ("`deyx'"=="") {
						if (!`fd') local rti "Predictive margins"
						else local rti "Contrasts of predictive margins"
						if ("`varlist'"!="") mata: __s_rlist = st_global("r(margins)")
					}
					else {
						local rti "Average marginal effects"
						mata: __s_rdyx = st_global("r(xvars)")
					}
					local rpl `r(predict_label)'
					local rex `r(expression)'
				}
				if ("`nonotes'"=="") local drv `r(derivatives)'

				cap confirm e `r(atstats1)'
				if (!_rc) {
					mata: __s_atm = st_matrix("r(at)")
					mata: __s_slb = tokens(st_global("r(atstats1)"))
				}

				if ("`over'"!="") {
					if (`xno'>1) {
						di as err "only one {it:xofivar} value allowed with {bf:over}"
						exit(198)
					}
					mata: ovarg(__s_acat,__s_aclv,__s_acname,__s_rac,"`esamp'")
				}
 
				if (`fd') {
					di as err "without {bf:nunit}, {bf:firstdiff} not allowed with {bf:xofinterest}"
					exit(198)
				}

				if ("`at'"!="") {
					mata: __s_selv = .
					mata: at_xac(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac,__s_selv) // most args are for [r]oarg()
					if ("`nolegend'"=="") {
						mata: __s_clb = st_matrixcolstripe("r(at)")
						mata: __s_atlv=__s_atname=__s_ats=__s_rat = .
						mata: atarg(__s_atlv,__s_atm,__s_atname,__s_ats,__s_atoz,__s_fvar,__s_rat,__s_selv)
					}
					else local nat = rowsof(r(at))/(`xno'*`nac')
				}
				else if (`nac'>1 & "`over'"=="") {
					if (`xno'>1) {
						di as err "only one {it:xofivar} value allowed with {bf:across}"
						exit(198)
					}
					mata: roarg(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac)
				}
				local nq = `nac'*`xno'
				if (`nq'==2) local rng 0

				// # of requested results (two elements)
				if (`nq') local nrq `nq'
				else local nrq 1
				if (`nac'==1) {
					local nres1 `nrq'
					local nres2 `nat'
				}
				else {
					local nres1 `nac'
					local nres2 `nat'
				}
				if ("`difference'"=="") {
					if (`rng' | `nres1'==2) local nres = `nres1'*`nres2'
					else local nres = (2*(`nres1'*`=`nres1'-1'/2))*`nres2'
				}	
				else {
					if (`rng' | `nres1'==2) local nres = `nres1'*`nres2'+`nres2'
					else local nres = (3*(`nres1'*`=`nres1'-1'/2))*`nres2'
				}
				if (`nres'>1000) {
					di as err "too many individual results requested; the limit is 1,000."
					exit(198)
				}
				else if (`nres'>100 & "`many'"=="") {
					di as err "more than 100 individual results requested; option {bf:many} required."
					exit(198)
				}

				local mc = ((`nq'>2 & !`rng') | `nac'>2 | `nat'>1)

				mata: __s_id=__s_ind=__s_rb=__s_rV = .
				mata: dsxarg(__s_id,__s_ind,__s_ord,__s_rb,__s_rV)
				
				if (!`rng') mata: sdi1(__s_df,__s_id,__s_ind,__s_lab,__s_matn,__s_rlab,__s_rb,__s_rV,__s_typ)
				else mata: sdi2(__s_df,__s_id,__s_ind,__s_lab,__s_matn,__s_rlab,__s_rb,__s_rV,__s_typ)	
			}
			else {
				if ("`xofinterest'"=="") {
					di as err "{bf:nunit} requires {bf:xofinterest}"
					exit(198)
				}
				if ((`nac'>1 | "`over'"!="") & !`fd') {
					di as err "with {bf:nunit}, {bf:across} requires {bf:firstdiff}"
					exit(198)
				}

				cap confirm number `nunit'
				if (_rc) {
					if ustrlower("`nunit'")=="sd" | ustrlower("`nunit'")=="-sd" {
						if length("`nunit'")==2 local sd 1
						else local sd 2				
					}
					else {
						di as err "invalid {bf:nunit} option"
						exit(198)
					}
				}
				else local sd 3

				margins in `fsamp', at(`xofinterest') nose noatlegend
				tempname xval1 xval2
				mata: __s_rx=__s_xat=__s_xlev=__s_xname = .		
				mata: xvarg(__s_acoz,__s_atoz,__s_naso,__s_ovoz,__s_rx,__s_xat,__s_xlev,__s_xname,__s_xvoz,"`esamp'","`xval1'","`xval2'")

				margins `if' `in' `fw' `aw' `pw' `iw', at(`=`acro'' `at' `=`xval1'') at(`=`acro'' `at' `=`xval2'') level(`level') over(`over') `options'
				local df = (el(r(table),7,1)!=.)
				mata: __s_df = st_matrix("r(table)")[7,1]
				mata: __s_rN = st_numscalar("r(N)")
				if ("`nolegend'"=="") {
					local len `=strlen("r(N)")'
					if (`len'>3) local rno `=string(`r(N)',"%-"+"`=`len'+int(`len'/3)+2'"+".0gc")'
					else local rno `=string(r(N))'
					if (!`fd') local rti = "Predictive margins"
					else local rti = "Contrasts of predictive margins"
					local rpl = "`r(predict_label)'"
					local rex = "`r(expression)'"
				}
				mata: __s_atm = st_matrix("r(at)")[1::rows(st_matrix("r(at)"))/2,]
				mata: __s_slb = tokens(st_global("r(atstats1)"))

				if ("`over'"!="") {
					if (`xno'>1) {
						di as err "only one {it:xofivar} value allowed with {bf:across}"
						exit(198)
					}
					mata: ovarg(__s_acat,__s_aclv,__s_acname,__s_rac,"`esamp'")
				}
			
				if ("`at'"!="") {
					mata: __s_selv = .
					mata: at_xac(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac,__s_selv) // most args are for [r]oarg()
					if ("`nolegend'"=="") {
						mata: __s_clb = st_matrixcolstripe("r(at)")
						mata: __s_atlv=__s_atname=__s_ats=__s_rat = .
						mata: atarg(__s_atlv,__s_atm,__s_atname,__s_ats,__s_atoz,__s_fvar,__s_rat,__s_selv)
					}
					else local nat = rowsof(r(at))/(2*`xno'*`nac')
				}
				else if (`nac'>1 & "`over'"=="") {
					if (`xno'>1) {
						di as err "only one {it:xofivar} value allowed with {bf:across}"
						exit(198)
					}
					mata: roarg(__s_acat,__s_acind,__s_aclv,__s_acname,__s_atm,__s_naso,__s_ord,__s_rac)
				}

				local nx = `xno'*`nac'	
				if ("`difference'"=="") { // # of requested results
					if (!`fd') local nres = 2*`nx'*`nat'
					else if (`rng' | `nx'==2) local nres = `nx'*`nat'
					else local nres = 2*(`nx'*`=`nx'-1'/2)*`nat'
				}	
				else {
					if (!`fd') local nres = (2*`nx'*`nat')+(`nx'*`nat')
					else if (`rng' | `nx'==2) local nres = `nx'*`nat'+`nat'
					else local nres = 3*(`nx'*`=`nx'-1'/2)*`nat'
				}
				if (`nres'>1000) {
					di as err "too many individual results requested; the limit is 1,000."
					exit(198)
				}
				else if (`nres'>100 & "`many'"=="") {
					noi di as text "more than 100 individual results requested; option {bf:many} required."
					exit(198)
				}

				local mc 0
				local lab "_xofi"
				if (`nat'==1) {
					if (!`fd' & `xno'>1) {
						local lab "Nº: _xofi "
						local mc 1 // multiple comparisons
					}
					else if (`fd') {
						if (`nac'==2 | (`nac'==2 & `rng')) {
							local lab "_acr#_xofi"
						}
						else if (`nac'>2) {
							local lab "Nº: _acr#_xofi "
							local mc 1
						}
						else {
							if (`xno'==2 | `rng') local lab "_xofi"
							else {
								local lab "Nº: _xofi "
								local mc 1
							}
						}
					}
				}
				else {
					if (!`fd') {
						local lab "Nº: _at#_xofi "
						local mc 1
					}				
					else {
						if (`nac'>1) {
							local lab "Nº: _at#_acr#_xofi "
							local mc 1
						}
						else {
							local lab "Nº: _at#_xofi "
							local mc 1
						}
					}
				}

				if (!`fd') mata: sdi3(__s_df,__s_lab,__s_matn,__s_ord,__s_rlab,__s_typ)
				else if (!`rng' | `nac'==2 | `xno'==2) mata: sdi4(__s_df,__s_lab,__s_matn,__s_ord,__s_rlab,__s_typ)
				else mata: sdi5(__s_df,__s_lab,__s_matn,__s_ord,__s_rlab,__s_typ)					
			}

			// Results presentation
			tempvar	acat aclv acname atlv atname ats ceq hi_sdi levtype lo_sdi ///
					quant ra rb rc statis std_err teq xat xeq xlev xname
			tempname lbl
			mata: st_vlmodify("`lbl'", __s_rlab, __s_lab)
			mata: st_local("nlab",strofreal(rows(__s_rlab)))

			noi di as text _newline "{bf:SDI} Results (pairwise comparisons)"
			
			mata: __s_gm = __s_matn[,(1,2,6,7)] // submatrix for 'getmata'
			getmata `quant'=__s_rlab `levtype'=__s_typ (`statis' `std_err' `lo_sdi' `hi_sdi')=__s_gm, force replace
			if ("`nolegend'"=="") {
				if ("`at'"!="") {
					getmata `ats'=__s_ats `atname'=__s_atname `atlv'=__s_atlv, force
					gen `teq' = "=" in 1/`ato'
					replace `teq' = ":" if substr(`atlv',1,1)=="("
				
					local lat =strlen("`nat'")+6
					format `ats' %-`lat's
					local len = strlen(`atname'[1])
					if (`len'>15) local len=15
					format `atname' %-`len's
					local fmt: format `atlv'
					local fmt: subinstr local fmt "%" "%-"
					format `atlv' `fmt'
					format `teq' %-1s
				}

				if (`nac'>1) {
					getmata `acat'=__s_acat `acname'=__s_acname `aclv'=__s_aclv, force
					gen `ceq' = "=" in 1/`nac'

					local lat = strlen("`nac'")+10
					format `acat' %-`lat's
					local len = strlen(`acname'[1])
					if (`len'>15) local len = 15
					format `acname' %-`len's
					local fmt: format `aclv'
					local fmt: subinstr local fmt "%" "%-"
					format `aclv' `fmt'
					format `ceq' %-1s
				}

				if ("`xofinterest'"!="") {
					getmata `xat'=__s_xat `xname'=__s_xname `xlev'=__s_xlev, force replace
					gen `xeq' = "=" in 1/`xno'
					replace `xeq' = ":" if substr(`xlev',1,1)=="("

					local lat = strlen("`xno'")+8
					format `xat' %-`lat's
					local len = strlen(`xname'[1]) 
					if (`len'>15) local len = 15 
					format `xname' %-`len's
					local fmt: format `xlev'
					local fmt: subinstr local fmt "%" "%-"
					format `xlev' `fmt'
					format `xeq' %-1s
				}

				gen `ra' = ""
				gen `rb' = `ra'
				gen `rc' = `ra'
				format `ra' %-15s
				format `rb' %-31s
				format `rc' %-1s
				
				replace `ra' = "Expression" in 1
				replace `rb' = "`rpl'"+", "+"`rex'" in 1
				replace `ra' = "Statistic" in 2
				replace `rb' = "`rti'" in 2
				replace `ra' = "Standard errors" in 3
				replace `rb' = "Delta-method" in 3
				if ("`mc'"!="1") {
					replace `ra' = "Number of obs" in 4
					replace `rb' = "`rno'" in 4
					replace	`rc' = ":" in 1/3
					replace `rc' = "=" in 4
					local nra 4
				}
				else {
					replace `ra' = "Nº" in 4
					replace `rb' = "Comparison reference number" in 4
					replace `ra' = "Number of obs" in 5
					replace `rb' = "`rno'" in 5				
					replace `rc' = ":" in 1/4
					replace `rc' = "=" in 5
					local nra 5
				}
				noi di as text ""
				noi flist `ra' `rc' `rb' in 1/`nra', noobs noh clean

				if ("`at'"!="") {
					noi di as text ""
					noi flist `ats' `atname' `teq' `atlv' in 1/`ato', noobs noh clean
				}
				if (`nac'>1) {
					noi di as text ""
					if ("`over'"=="") noi di as text _skip(4) "{bf:acr[oss]}"
					if ("`over'"!="") noi di as text _skip(4) "{bf:acr[oss](over)}"
					noi flist `acat' `acname' `ceq' `aclv' in 1/`nac', noobs noh clean
				}
				if ("`xofinterest'"!="") {
					noi di as text "" _n _skip(4) "{bf:xofi[nterest]}"
					noi flist `xat' `xname' `xeq' `xlev' in 1/`xno', noobs noh clean
					if ("`nunit'"!="") {
						if ("`sd'"=="3") noi di as text _newline _skip(4) "{it:n}{bf:-unit   =   `unit'}"
						else noi di as text _newline _skip(4) "{it:n}{bf:-unit   =   `unit' {bf:(std. dev.)}}"
					}
				}
			}

			label var `quant' "`=char(13)'"
			label values `quant' `lbl'
			label var `statis' "Statistic"
			label var `std_err' "Std. Err."
			label var `lo_sdi' "[Interval"
			label var `hi_sdi' "Bounds]"
			label var `levtype' "(%) Type"

			noi tabdis `quant' in 1/`nlab', c(`statis' `std_err' `lo_sdi' `hi_sdi' `levtype') left cellwidth(13)
			if ("`nonotes'"=="") noi di as text "Note: SDIs indicate significance of difference from `mvalue'."

			mata: st_matrix("r(sdi)",__s_matn[3::rows(__s_matn),])
			mata: __s_lab = __s_lab[3::rows(__s_lab),]
			if (`mc') mata: __s_lab = strtrim(substr(__s_lab,1,strpos(__s_lab,":"):-1):+"._no")
			else mata: __s_lab = J(rows(__s_lab),1,"1._no")
			mata: st_matrixrowstripe("r(sdi)",(J(rows(__s_lab),1,""),__s_lab))
			if (!`df') {
				if (`cens') mata: st_matrixcolstripe("r(sdi)",(J(8,1,""),("b"\"se"\"z"\"pvalue"\"c(crit)"\"ll"\"ul"\"level")))
				else mata: st_matrixcolstripe("r(sdi)",(J(8,1,""),("b"\"se"\"z"\"pvalue"\"crit"\"ll"\"ul"\"level")))
			}
			else {
				if (`cens') mata: st_matrixcolstripe("r(sdi)",(J(8,1,""),("b"\"se"\"t"\"pvalue"\"c(crit)"\"ll"\"ul"\"level")))
				else mata: st_matrixcolstripe("r(sdi)",(J(8,1,""),("b"\"se"\"t"\"pvalue"\"crit"\"ll"\"ul"\"level")))
			}
		}
	}

	if (_rc) {
		local rc = _rc
		cap if (`N'>_N) drop in `=`N'+1'/`=_N'
		cap mata: mata drop __s_*
		exit(`rc')
	}
	else {
		cap if (`N'>_N) drop in `=`N'+1'/`=_N'
		cap mata: mata drop __s_*
	}

end
