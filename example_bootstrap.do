use data/df_hom.dta, clear

egen unique_id = group(state unit)

capture program drop did2s_est

program did2s_est, rclass
	version 13.0
	regress dep_var i.new_id i.year if treat == 0
	tempvar dep_var_resid
	predict `dep_var_resid', residuals
	
	regress `dep_var_resid' ib0.treat

end

xtset unique_id year
sort unique_id year
bootstrap, cluster(state) idcluster(new_id) group(unique_id) reps(100): did2s_est
