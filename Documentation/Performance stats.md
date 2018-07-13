20/6-18:

get_negL for 10 000 GDR1 stars when calc_K has a for-loop:

	3.2 s 	for n = 8x8x8, dv = 50x50x50, vmin = (-200,-200,200)
	
	23.1 s	for n = 50x50x50, dv = 8x8x8, vmin = (-200,-200,-200)
	
get_negL for 10 000 GDR1 stars when calc_K just use arrays:

	2.6 s	for n = 8x8x8, dv = 50x50x50, vmin = (-200,-200,-200)
	
	16.8 s	for n = 50x50x50, dv = 8x8x8, vmin = (-200,-200,-200)
	
13/07/18:

	get_negL for N = 2000, n = 10x10x10, dv = (40,40,40), vmin = -200
	
		33 function calls in 0.013 seconds

		   Ordered by: standard name

	   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
			1    0.002    0.002    0.013    0.013 <string>:1(<module>)
			1    0.000    0.000    0.000    0.000 Deproject_v0.py:140(sec_der)
			1    0.008    0.008    0.011    0.011 Deproject_v0.py:216(get_negL)
			5    0.000    0.000    0.002    0.000 _methods.py:31(_sum)
			1    0.000    0.000    0.000    0.000 fromnumeric.py:138(reshape)
			5    0.000    0.000    0.003    0.001 fromnumeric.py:1730(sum)
			1    0.000    0.000    0.000    0.000 fromnumeric.py:55(_wrapfunc)
			1    0.000    0.000    0.013    0.013 {built-in method builtins.exec}
			1    0.000    0.000    0.000    0.000 {built-in method builtins.getattr}
			5    0.000    0.000    0.000    0.000 {built-in method builtins.isinstance}
			1    0.000    0.000    0.000    0.000 {built-in method builtins.len}
			1    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.array}
			1    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.zeros}
			1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
			5    0.002    0.000    0.002    0.000 {method 'reduce' of 'numpy.ufunc' objects}
			2    0.000    0.000    0.000    0.000 {method 'reshape' of 'numpy.ndarray' objects}
			
	get_negL for N = 2000, n = 10x10x10, dv = (40,40,40), vmin = -200
	
120718 results:

	Best run so far

	N = 1500
	v0 = 0,0,0
	disp = 50,50,50
	n = 25x25x25
	dv = 16x16x16
	v0_guess = 2,1,-3
	disp_guess = 30,30,30
	
	Sample mean ===> [-0.25048585  0.65953818 -2.12676714]
	Sample velocity dispersion ===> [ 50.45938417  50.68601034  49.58213394]
	Computed mean ===> [-0.07797956  0.62970625 -2.71046899]
	Computed velocity dispersion ===> [ 48.46015048  47.12892852  46.47691713]
	
max_L runtimes:

	cProfile stats for max_L, N = 2000, n = 10x10x10, dv = (40,40,40), gtol = 1e-4

	2631881 function calls (2631869 primitive calls) in 58.324 seconds

	   Ordered by: standard name

	   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
			1    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:997(_handle_fromlist)
			1    0.002    0.002   58.324   58.324 <string>:1(<module>)
			9    0.000    0.000    0.000    0.000 <string>:12(__new__)
			1    0.000    0.000    0.000    0.000 <string>:2(<module>)
			2    0.000    0.000    0.000    0.000 <string>:2(_parse_args)
			1    0.001    0.001    0.008    0.008 Deproject_v0.py:107(calc_sigma2)
		 1418    0.157    0.000    0.218    0.000 Deproject_v0.py:140(sec_der)
			1    0.000    0.000    0.006    0.006 Deproject_v0.py:179(phi_guess)
		  715    7.431    0.010    9.787    0.014 Deproject_v0.py:216(get_negL)
		  703   14.524    0.021   41.927    0.060 Deproject_v0.py:245(get_grad_negL)
			1    0.032    0.032   58.322   58.322 Deproject_v0.py:277(max_L)
		 2000    0.208    0.000    0.538    0.000 Deproject_v0.py:44(calc_K)
			3    0.000    0.000    0.000    0.000 _continuous_distns.py:137(_logpdf)
			3    0.000    0.000    0.000    0.000 _continuous_distns.py:80(_norm_logpdf)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:1484(__init__)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:1549(_updated_ctor_param)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:1663(logpdf)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:429(__init__)
			9    0.000    0.000    0.000    0.000 _distn_infrastructure.py:43(instancemethod)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:452(logpdf)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:524(argsreduce)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:549(<listcomp>)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:591(__init__)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:626(_construct_argparser)
			3    0.000    0.000    0.003    0.001 _distn_infrastructure.py:706(_construct_doc)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:714(<genexpr>)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:754(freeze)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:771(__call__)
			6    0.000    0.000    0.000    0.000 _distn_infrastructure.py:863(_argcheck)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:875(_support_mask)
		  134    0.000    0.000    0.003    0.000 _methods.py:25(_amax)
		 6388    0.009    0.000    5.598    0.001 _methods.py:31(_sum)
			3    0.000    0.000    0.000    0.000 _methods.py:37(_any)
			3    0.000    0.000    0.000    0.000 _methods.py:43(_count_reduce_items)
			3    0.000    0.000    0.001    0.000 _methods.py:53(_mean)
			3    0.000    0.000    0.000    0.000 _util.py:173(check_random_state)
			9    0.000    0.000    0.001    0.000 _util.py:269(getargspec_no_self)
			9    0.000    0.000    0.000    0.000 _util.py:292(<listcomp>)
			9    0.000    0.000    0.000    0.000 _util.py:296(<listcomp>)
			9    0.000    0.000    0.000    0.000 _util.py:301(<listcomp>)
			9    0.000    0.000    0.000    0.000 _util.py:306(<listcomp>)
			6    0.001    0.000    0.002    0.000 doccer.py:12(docformat)
			6    0.000    0.000    0.001    0.000 doccer.py:128(indentcount_lines)
		 1418    0.005    0.000    0.023    0.000 fromnumeric.py:138(reshape)
		  716    0.003    0.000    0.008    0.000 fromnumeric.py:1380(ravel)
		 2006    0.002    0.000    0.007    0.000 fromnumeric.py:1487(nonzero)
			3    0.000    0.000    0.000    0.000 fromnumeric.py:1565(shape)
		 6388    0.056    0.000    5.661    0.001 fromnumeric.py:1730(sum)
			3    0.000    0.000    0.000    0.000 fromnumeric.py:1886(any)
		  134    0.001    0.000    0.004    0.000 fromnumeric.py:2174(amax)
			3    0.000    0.000    0.001    0.000 fromnumeric.py:2806(mean)
		 3430    0.006    0.000    0.023    0.000 fromnumeric.py:55(_wrapfunc)
			6    0.000    0.000    0.000    0.000 fromnumeric.py:70(take)
			3    0.000    0.000    0.000    0.000 function_base.py:213(iterable)
			6    0.000    0.000    0.000    0.000 function_base.py:2283(extract)
			3    0.000    0.000    0.000    0.000 function_base.py:2334(place)
		   12    0.000    0.000    0.000    0.000 function_base.py:2679(__init__)
			1    0.000    0.000    0.000    0.000 function_base.py:4554(meshgrid)
			1    0.000    0.000    0.000    0.000 function_base.py:4671(<listcomp>)
			1    0.000    0.000    0.000    0.000 function_base.py:4684(<listcomp>)
		   18    0.000    0.000    0.000    0.000 inspect.py:159(isfunction)
			9    0.000    0.000    0.000    0.000 inspect.py:1777(_signature_bound_method)
			9    0.000    0.000    0.000    0.000 inspect.py:2092(_signature_from_function)
		 18/9    0.000    0.000    0.001    0.000 inspect.py:2173(_signature_from_callable)
		   15    0.000    0.000    0.000    0.000 inspect.py:2431(__init__)
		   27    0.000    0.000    0.000    0.000 inspect.py:2480(name)
		   12    0.000    0.000    0.000    0.000 inspect.py:2484(default)
		   48    0.000    0.000    0.000    0.000 inspect.py:2492(kind)
		   18    0.000    0.000    0.000    0.000 inspect.py:2710(__init__)
		   24    0.000    0.000    0.000    0.000 inspect.py:2755(<genexpr>)
			9    0.000    0.000    0.001    0.000 inspect.py:2779(from_callable)
		   45    0.000    0.000    0.000    0.000 inspect.py:2785(parameters)
			9    0.000    0.000    0.000    0.000 inspect.py:2793(replace)
			9    0.000    0.000    0.001    0.000 inspect.py:3031(signature)
			9    0.000    0.000    0.000    0.000 inspect.py:485(unwrap)
			9    0.000    0.000    0.000    0.000 inspect.py:505(_is_wrapper)
		   11    0.000    0.000    0.000    0.000 iostream.py:180(schedule)
		   10    0.000    0.000    0.000    0.000 iostream.py:284(_is_master_process)
		   10    0.000    0.000    0.000    0.000 iostream.py:297(_schedule_flush)
		   10    0.000    0.000    0.000    0.000 iostream.py:342(write)
		   11    0.000    0.000    0.000    0.000 iostream.py:87(_event_pipe)
			1    0.000    0.000    0.000    0.000 linalg.py:101(get_linalg_error_extobj)
			1    0.000    0.000    0.000    0.000 linalg.py:106(_makearray)
			3    0.000    0.000    0.000    0.000 linalg.py:111(isComplexType)
			1    0.000    0.000    0.000    0.000 linalg.py:124(_realType)
			1    0.000    0.000    0.000    0.000 linalg.py:139(_commonType)
			1    0.000    0.000    0.000    0.000 linalg.py:198(_assertRankAtLeast2)
			1    0.000    0.000    0.000    0.000 linalg.py:2014(norm)
			1    0.000    0.000    0.000    0.000 linalg.py:209(_assertNdSquareness)
			1    0.000    0.000    0.000    0.000 linalg.py:449(inv)
		  134    0.018    0.000   57.473    0.429 linesearch.py:106(scalar_search_wolfe1)
			1    0.000    0.000    0.167    0.167 linesearch.py:194(line_search_wolfe2)
		   12    0.000    0.000    0.165    0.014 linesearch.py:257(phi)
			1    0.000    0.000    0.167    0.167 linesearch.py:296(scalar_search_wolfe2)
		  134    0.001    0.000   57.475    0.429 linesearch.py:34(line_search_wolfe1)
		   10    0.001    0.000    0.001    0.000 linesearch.py:419(_cubicmin)
			1    0.000    0.000    0.000    0.000 linesearch.py:453(_quadmin)
			1    0.000    0.000    0.153    0.153 linesearch.py:474(_zoom)
		  702    0.020    0.000   11.125    0.016 linesearch.py:85(phi)
		  702    0.025    0.000   46.329    0.066 linesearch.py:89(derphi)
			1    0.000    0.000    0.000    0.000 numeric.py:2365(identity)
		   22    0.000    0.000    0.000    0.000 numeric.py:2667(seterr)
		   22    0.000    0.000    0.000    0.000 numeric.py:2767(geterr)
		   11    0.000    0.000    0.000    0.000 numeric.py:3060(__init__)
		   11    0.000    0.000    0.000    0.000 numeric.py:3064(__enter__)
		   11    0.000    0.000    0.000    0.000 numeric.py:3069(__exit__)
		   25    0.000    0.000    0.000    0.000 numeric.py:463(asarray)
	   823331    0.188    0.000    0.449    0.000 numeric.py:534(asanyarray)
			3    0.000    0.000    0.000    0.000 numerictypes.py:1015(<listcomp>)
			3    0.000    0.000    0.000    0.000 numerictypes.py:1016(<listcomp>)
			6    0.000    0.000    0.000    0.000 numerictypes.py:942(_can_coerce_all)
		   42    0.000    0.000    0.000    0.000 numerictypes.py:951(<listcomp>)
			3    0.000    0.000    0.000    0.000 numerictypes.py:964(find_common_type)
			1    0.000    0.000   57.737   57.737 optimize.py:1020(fmin_cg)
			1    0.009    0.009   57.737   57.737 optimize.py:1191(_minimize_cg)
			1    0.000    0.000    0.000    0.000 optimize.py:137(_check_unknown_options)
		  134    0.001    0.000    0.005    0.000 optimize.py:155(vecnorm)
			2    0.000    0.000    0.000    0.000 optimize.py:285(wrap_function)
		 1418    5.924    0.004   57.638    0.041 optimize.py:290(function_wrapper)
		  134    0.001    0.000   57.643    0.430 optimize.py:750(_line_search_wolfe12)
			3    0.000    0.000    0.000    0.000 shape_base.py:11(atleast_1d)
		 8704    0.208    0.000   24.124    0.003 shape_base.py:296(stack)
		 8704    0.250    0.000    0.697    0.000 shape_base.py:348(<listcomp>)
	   831304    0.145    0.000    0.145    0.000 shape_base.py:352(<genexpr>)
		 8704    0.263    0.000    0.263    0.000 shape_base.py:360(<listcomp>)
			3    0.000    0.000    0.000    0.000 stride_tricks.py:115(_broadcast_to)
		   12    0.000    0.000    0.000    0.000 stride_tricks.py:120(<genexpr>)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:176(_broadcast_shape)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:195(broadcast_arrays)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:247(<listcomp>)
			3    0.000    0.000    0.000    0.000 stride_tricks.py:25(_maybe_view_as_subclass)
			2    0.000    0.000    0.000    0.000 stride_tricks.py:251(<genexpr>)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:257(<listcomp>)
		   11    0.000    0.000    0.000    0.000 threading.py:1062(_wait_for_tstate_lock)
		   11    0.000    0.000    0.000    0.000 threading.py:1104(is_alive)
		   11    0.000    0.000    0.000    0.000 threading.py:506(is_set)
			1    0.000    0.000    0.000    0.000 twodim_base.py:139(eye)
			1    0.000    0.000    0.000    0.000 warnings.py:143(simplefilter)
			1    0.000    0.000    0.000    0.000 warnings.py:159(_add_filter)
			1    0.000    0.000    0.000    0.000 warnings.py:428(__init__)
			1    0.000    0.000    0.000    0.000 warnings.py:449(__enter__)
			1    0.000    0.000    0.000    0.000 warnings.py:468(__exit__)
			9    0.000    0.000    0.000    0.000 {built-in method __new__ of type object at 0x000000006BD9B3F0}
			3    0.000    0.000    0.000    0.000 {built-in method _warnings._filters_mutated}
			1    0.000    0.000    0.000    0.000 {built-in method _warnings.warn}
			1    0.000    0.000    0.000    0.000 {built-in method builtins.all}
			3    0.000    0.000    0.000    0.000 {built-in method builtins.any}
		   18    0.000    0.000    0.000    0.000 {built-in method builtins.callable}
		  4/1    0.001    0.000   58.324   58.324 {built-in method builtins.exec}
		 3431    0.002    0.000    0.002    0.000 {built-in method builtins.getattr}
		   14    0.000    0.000    0.000    0.000 {built-in method builtins.hasattr}
			9    0.000    0.000    0.000    0.000 {built-in method builtins.id}
		 7334    0.008    0.000    0.008    0.000 {built-in method builtins.isinstance}
		   11    0.000    0.000    0.000    0.000 {built-in method builtins.issubclass}
			3    0.000    0.000    0.000    0.000 {built-in method builtins.iter}
		19478    0.006    0.000    0.006    0.000 {built-in method builtins.len}
		 8134    0.018    0.000    0.018    0.000 {built-in method builtins.max}
		 8430    0.015    0.000    0.015    0.000 {built-in method builtins.min}
			5    0.000    0.000    0.000    0.000 {built-in method builtins.print}
			9    0.000    0.000    0.000    0.000 {built-in method builtins.setattr}
		   10    0.000    0.000    0.000    0.000 {built-in method nt.getpid}
		   11    0.000    0.000    0.000    0.000 {built-in method nt.urandom}
			3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray._insert}
		 6003    0.011    0.000    0.011    0.000 {built-in method numpy.core.multiarray.arange}
	   824781    0.285    0.000    0.285    0.000 {built-in method numpy.core.multiarray.array}
		10704   22.811    0.002   22.811    0.002 {built-in method numpy.core.multiarray.concatenate}
		 1118    0.016    0.000    0.016    0.000 {built-in method numpy.core.multiarray.dot}
		   13    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.empty}
		 8704    0.004    0.000    0.004    0.000 {built-in method numpy.core.multiarray.normalize_axis_index}
			3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.putmask}
		 3688    0.017    0.000    0.017    0.000 {built-in method numpy.core.multiarray.zeros}
		   44    0.000    0.000    0.000    0.000 {built-in method numpy.core.umath.geterrobj}
		   22    0.000    0.000    0.000    0.000 {built-in method numpy.core.umath.seterrobj}
			1    0.000    0.000    0.000    0.000 {method '__array_prepare__' of 'numpy.ndarray' objects}
		   11    0.000    0.000    0.000    0.000 {method 'acquire' of '_thread.lock' objects}
			3    0.000    0.000    0.000    0.000 {method 'any' of 'numpy.ndarray' objects}
		 1612    0.000    0.000    0.000    0.000 {method 'append' of 'list' objects}
		 2001    0.003    0.000    0.003    0.000 {method 'astype' of 'numpy.ndarray' objects}
			6    0.000    0.000    0.000    0.000 {method 'copy' of 'dict' objects}
			3    0.000    0.000    0.000    0.000 {method 'copy' of 'numpy.ndarray' objects}
			1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
		  192    0.000    0.000    0.000    0.000 {method 'expandtabs' of 'str' objects}
			3    0.000    0.000    0.000    0.000 {method 'fill' of 'numpy.ndarray' objects}
		   11    0.000    0.000    0.000    0.000 {method 'flatten' of 'numpy.ndarray' objects}
		   28    0.000    0.000    0.000    0.000 {method 'get' of 'dict' objects}
			1    0.000    0.000    0.000    0.000 {method 'insert' of 'list' objects}
		   15    0.000    0.000    0.000    0.000 {method 'isidentifier' of 'str' objects}
			6    0.000    0.000    0.000    0.000 {method 'items' of 'dict' objects}
		  165    0.000    0.000    0.000    0.000 {method 'join' of 'str' objects}
		  393    0.000    0.000    0.000    0.000 {method 'lstrip' of 'str' objects}
		 2006    0.002    0.000    0.002    0.000 {method 'nonzero' of 'numpy.ndarray' objects}
		   15    0.000    0.000    0.000    0.000 {method 'pop' of 'dict' objects}
		  717    0.003    0.000    0.003    0.000 {method 'ravel' of 'numpy.ndarray' objects}
		 6528    5.593    0.001    5.593    0.001 {method 'reduce' of 'numpy.ufunc' objects}
			1    0.000    0.000    0.000    0.000 {method 'remove' of 'list' objects}
		   18    0.000    0.000    0.000    0.000 {method 'replace' of 'str' objects}
		 3546    0.027    0.000    0.027    0.000 {method 'reshape' of 'numpy.ndarray' objects}
		 2000    0.004    0.000    0.004    0.000 {method 'sort' of 'numpy.ndarray' objects}
		  192    0.000    0.000    0.000    0.000 {method 'splitlines' of 'str' objects}
			6    0.000    0.000    0.000    0.000 {method 'take' of 'numpy.ndarray' objects}
		   45    0.000    0.000    0.000    0.000 {method 'values' of 'mappingproxy' objects}
		   
	
After optimising by creating arrays and assigning values instead of stacking we got and using gtol = 5e-4:
	
	110907 function calls (110895 primitive calls) in 28.069 seconds

	Ordered by: standard name

	   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
			1    0.000    0.000    0.000    0.000 <frozen importlib._bootstrap>:997(_handle_fromlist)
			1    0.002    0.002   28.068   28.068 <string>:1(<module>)
			9    0.000    0.000    0.000    0.000 <string>:12(__new__)
			1    0.000    0.000    0.000    0.000 <string>:2(<module>)
			2    0.000    0.000    0.000    0.000 <string>:2(_parse_args)
			1    0.000    0.000    0.002    0.002 Deproject_v0.py:111(calc_sigma2)
		 1040    0.116    0.000    0.160    0.000 Deproject_v0.py:148(sec_der)
			1    0.000    0.000    0.006    0.006 Deproject_v0.py:187(phi_guess)
		  520    5.375    0.010    7.107    0.014 Deproject_v0.py:224(get_negL)
		  520   12.916    0.025   15.804    0.030 Deproject_v0.py:253(get_grad_negL)
			1    0.047    0.047   28.066   28.066 Deproject_v0.py:287(max_L)
		 2000    0.203    0.000    0.289    0.000 Deproject_v0.py:44(calc_K)
			3    0.000    0.000    0.000    0.000 _continuous_distns.py:145(_logpdf)
			3    0.000    0.000    0.000    0.000 _continuous_distns.py:85(_norm_logpdf)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:1487(__init__)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:1552(_updated_ctor_param)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:1666(logpdf)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:429(__init__)
			9    0.000    0.000    0.000    0.000 _distn_infrastructure.py:43(instancemethod)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:452(logpdf)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:524(argsreduce)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:549(<listcomp>)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:591(__init__)
			3    0.000    0.000    0.001    0.000 _distn_infrastructure.py:626(_construct_argparser)
			3    0.000    0.000    0.003    0.001 _distn_infrastructure.py:706(_construct_doc)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:714(<genexpr>)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:754(freeze)
			3    0.000    0.000    0.004    0.001 _distn_infrastructure.py:771(__call__)
			6    0.000    0.000    0.000    0.000 _distn_infrastructure.py:863(_argcheck)
			3    0.000    0.000    0.000    0.000 _distn_infrastructure.py:875(_support_mask)
		  113    0.000    0.000    0.003    0.000 _methods.py:25(_amax)
		 4681    0.006    0.000    4.362    0.001 _methods.py:31(_sum)
			3    0.000    0.000    0.000    0.000 _methods.py:37(_any)
			3    0.000    0.000    0.000    0.000 _methods.py:43(_count_reduce_items)
			3    0.000    0.000    0.000    0.000 _methods.py:53(_mean)
			3    0.000    0.000    0.000    0.000 _util.py:173(check_random_state)
			9    0.000    0.000    0.001    0.000 _util.py:269(getargspec_no_self)
			9    0.000    0.000    0.000    0.000 _util.py:292(<listcomp>)
			9    0.000    0.000    0.000    0.000 _util.py:296(<listcomp>)
			9    0.000    0.000    0.000    0.000 _util.py:301(<listcomp>)
			9    0.000    0.000    0.000    0.000 _util.py:306(<listcomp>)
			6    0.001    0.000    0.002    0.000 doccer.py:12(docformat)
			6    0.000    0.000    0.001    0.000 doccer.py:172(indentcount_lines)
		   15    0.000    0.000    0.000    0.000 enum.py:265(__call__)
		   15    0.000    0.000    0.000    0.000 enum.py:515(__new__)
		  533    0.002    0.000    0.006    0.000 fromnumeric.py:1427(ravel)
		 2006    0.002    0.000    0.007    0.000 fromnumeric.py:1534(nonzero)
			3    0.000    0.000    0.000    0.000 fromnumeric.py:1613(shape)
		 1040    0.004    0.000    0.016    0.000 fromnumeric.py:163(reshape)
		 4681    0.041    0.000    4.409    0.001 fromnumeric.py:1778(sum)
			3    0.000    0.000    0.000    0.000 fromnumeric.py:1934(any)
		  113    0.001    0.000    0.004    0.000 fromnumeric.py:2222(amax)
			3    0.000    0.000    0.000    0.000 fromnumeric.py:2854(mean)
		 3052    0.005    0.000    0.018    0.000 fromnumeric.py:50(_wrapfunc)
			6    0.000    0.000    0.000    0.000 fromnumeric.py:65(take)
			3    0.000    0.000    0.000    0.000 function_base.py:213(iterable)
			6    0.000    0.000    0.000    0.000 function_base.py:2299(extract)
			3    0.000    0.000    0.000    0.000 function_base.py:2350(place)
		   12    0.000    0.000    0.000    0.000 function_base.py:2695(__init__)
			1    0.000    0.000    0.000    0.000 function_base.py:4568(meshgrid)
			1    0.000    0.000    0.000    0.000 function_base.py:4685(<listcomp>)
			1    0.000    0.000    0.000    0.000 function_base.py:4698(<listcomp>)
		   18    0.000    0.000    0.000    0.000 inspect.py:159(isfunction)
			9    0.000    0.000    0.000    0.000 inspect.py:1780(_signature_bound_method)
			9    0.000    0.000    0.000    0.000 inspect.py:2095(_signature_from_function)
		 18/9    0.000    0.000    0.001    0.000 inspect.py:2176(_signature_from_callable)
		   15    0.000    0.000    0.000    0.000 inspect.py:2445(__init__)
		   27    0.000    0.000    0.000    0.000 inspect.py:2495(name)
		   12    0.000    0.000    0.000    0.000 inspect.py:2499(default)
		   48    0.000    0.000    0.000    0.000 inspect.py:2507(kind)
		   18    0.000    0.000    0.000    0.000 inspect.py:2725(__init__)
		   24    0.000    0.000    0.000    0.000 inspect.py:2774(<genexpr>)
			9    0.000    0.000    0.001    0.000 inspect.py:2798(from_callable)
		   45    0.000    0.000    0.000    0.000 inspect.py:2804(parameters)
			9    0.000    0.000    0.000    0.000 inspect.py:2812(replace)
			9    0.000    0.000    0.001    0.000 inspect.py:3050(signature)
			9    0.000    0.000    0.000    0.000 inspect.py:485(unwrap)
			9    0.000    0.000    0.000    0.000 inspect.py:505(_is_wrapper)
		   11    0.000    0.000    0.002    0.000 iostream.py:195(schedule)
		   10    0.000    0.000    0.000    0.000 iostream.py:300(_is_master_process)
		   10    0.000    0.000    0.000    0.000 iostream.py:313(_schedule_flush)
		   10    0.000    0.000    0.002    0.000 iostream.py:366(write)
		   11    0.000    0.000    0.000    0.000 iostream.py:93(_event_pipe)
			1    0.000    0.000    0.000    0.000 linalg.py:100(get_linalg_error_extobj)
			1    0.000    0.000    0.000    0.000 linalg.py:105(_makearray)
			3    0.000    0.000    0.000    0.000 linalg.py:110(isComplexType)
			1    0.000    0.000    0.000    0.000 linalg.py:123(_realType)
			1    0.000    0.000    0.000    0.000 linalg.py:138(_commonType)
			1    0.000    0.000    0.000    0.000 linalg.py:197(_assertRankAtLeast2)
			1    0.000    0.000    0.000    0.000 linalg.py:208(_assertNdSquareness)
			1    0.000    0.000    0.000    0.000 linalg.py:2103(norm)
			1    0.000    0.000    0.000    0.000 linalg.py:464(inv)
		  112    0.014    0.000   27.632    0.247 linesearch.py:106(scalar_search_wolfe1)
		  112    0.001    0.000   27.634    0.247 linesearch.py:34(line_search_wolfe1)
		  519    0.014    0.000    8.311    0.016 linesearch.py:85(phi)
		  519    0.018    0.000   19.307    0.037 linesearch.py:89(derphi)
		 2000    0.002    0.000    0.011    0.000 numeric.py:146(ones)
			1    0.000    0.000    0.000    0.000 numeric.py:2157(identity)
		   15    0.000    0.000    0.000    0.000 numeric.py:424(asarray)
		  548    0.001    0.000    0.002    0.000 numeric.py:495(asanyarray)
			3    0.000    0.000    0.000    0.000 numerictypes.py:1013(<listcomp>)
			3    0.000    0.000    0.000    0.000 numerictypes.py:1014(<listcomp>)
			6    0.000    0.000    0.000    0.000 numerictypes.py:939(_can_coerce_all)
		   42    0.000    0.000    0.000    0.000 numerictypes.py:948(<listcomp>)
			3    0.000    0.000    0.000    0.000 numerictypes.py:962(find_common_type)
			1    0.000    0.000   27.722   27.722 optimize.py:1057(fmin_cg)
			1    0.001    0.001   27.722   27.722 optimize.py:1228(_minimize_cg)
		  112    0.004    0.000    0.010    0.000 optimize.py:1286(polak_ribiere_powell_step)
		  112    0.001    0.000    0.011    0.000 optimize.py:1296(descent_condition)
			1    0.000    0.000    0.000    0.000 optimize.py:137(_check_unknown_options)
		  113    0.001    0.000    0.004    0.000 optimize.py:156(vecnorm)
			2    0.000    0.000    0.000    0.000 optimize.py:286(wrap_function)
		 1040    4.735    0.005   27.647    0.027 optimize.py:291(function_wrapper)
		  112    0.004    0.000   27.648    0.247 optimize.py:785(_line_search_wolfe12)
			3    0.000    0.000    0.000    0.000 shape_base.py:11(atleast_1d)
		   11    0.001    0.000    0.001    0.000 socket.py:333(send)
			3    0.000    0.000    0.000    0.000 stride_tricks.py:115(_broadcast_to)
		   12    0.000    0.000    0.000    0.000 stride_tricks.py:120(<genexpr>)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:176(_broadcast_shape)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:195(broadcast_arrays)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:247(<listcomp>)
			3    0.000    0.000    0.000    0.000 stride_tricks.py:25(_maybe_view_as_subclass)
			2    0.000    0.000    0.000    0.000 stride_tricks.py:251(<genexpr>)
			1    0.000    0.000    0.000    0.000 stride_tricks.py:257(<listcomp>)
		   11    0.000    0.000    0.000    0.000 threading.py:1062(_wait_for_tstate_lock)
		   11    0.000    0.000    0.000    0.000 threading.py:1104(is_alive)
		   11    0.000    0.000    0.000    0.000 threading.py:506(is_set)
			1    0.000    0.000    0.000    0.000 twodim_base.py:140(eye)
			9    0.000    0.000    0.000    0.000 {built-in method __new__ of type object at 0x000000006BD9C3F0}
			1    0.000    0.000    0.000    0.000 {built-in method builtins.all}
			3    0.000    0.000    0.000    0.000 {built-in method builtins.any}
		   18    0.000    0.000    0.000    0.000 {built-in method builtins.callable}
		  4/1    0.001    0.000   28.069   28.069 {built-in method builtins.exec}
		 3053    0.002    0.000    0.002    0.000 {built-in method builtins.getattr}
		   14    0.000    0.000    0.000    0.000 {built-in method builtins.hasattr}
			9    0.000    0.000    0.000    0.000 {built-in method builtins.id}
		 5421    0.006    0.000    0.006    0.000 {built-in method builtins.isinstance}
		   11    0.000    0.000    0.000    0.000 {built-in method builtins.issubclass}
			3    0.000    0.000    0.000    0.000 {built-in method builtins.iter}
		18225    0.004    0.000    0.004    0.000 {built-in method builtins.len}
		 8113    0.017    0.000    0.017    0.000 {built-in method builtins.max}
		 8413    0.014    0.000    0.014    0.000 {built-in method builtins.min}
			5    0.000    0.000    0.002    0.000 {built-in method builtins.print}
			9    0.000    0.000    0.000    0.000 {built-in method builtins.setattr}
		   10    0.000    0.000    0.000    0.000 {built-in method nt.getpid}
			3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray._insert}
		 6003    0.011    0.000    0.011    0.000 {built-in method numpy.core.multiarray.arange}
		 1610    0.018    0.000    0.018    0.000 {built-in method numpy.core.multiarray.array}
		 2000    0.005    0.000    0.005    0.000 {built-in method numpy.core.multiarray.concatenate}
		 2000    0.006    0.000    0.006    0.000 {built-in method numpy.core.multiarray.copyto}
		 1081    0.011    0.000    0.011    0.000 {built-in method numpy.core.multiarray.dot}
		 2003    0.002    0.000    0.002    0.000 {built-in method numpy.core.multiarray.empty}
			3    0.000    0.000    0.000    0.000 {built-in method numpy.core.multiarray.putmask}
		 9787    0.057    0.000    0.057    0.000 {built-in method numpy.core.multiarray.zeros}
			9    0.000    0.000    0.000    0.000 {built-in method sys.getrecursionlimit}
			1    0.000    0.000    0.000    0.000 {method '__array_prepare__' of 'numpy.ndarray' objects}
		   11    0.000    0.000    0.000    0.000 {method 'acquire' of '_thread.lock' objects}
			3    0.000    0.000    0.000    0.000 {method 'any' of 'numpy.ndarray' objects}
		   11    0.000    0.000    0.000    0.000 {method 'append' of 'collections.deque' objects}
		 1591    0.000    0.000    0.000    0.000 {method 'append' of 'list' objects}
		 2001    0.003    0.000    0.003    0.000 {method 'astype' of 'numpy.ndarray' objects}
			6    0.000    0.000    0.000    0.000 {method 'copy' of 'dict' objects}
			3    0.000    0.000    0.000    0.000 {method 'copy' of 'numpy.ndarray' objects}
			1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
		  192    0.000    0.000    0.000    0.000 {method 'expandtabs' of 'str' objects}
			3    0.000    0.000    0.000    0.000 {method 'fill' of 'numpy.ndarray' objects}
			1    0.000    0.000    0.000    0.000 {method 'flatten' of 'numpy.ndarray' objects}
		   28    0.000    0.000    0.000    0.000 {method 'get' of 'dict' objects}
		   15    0.000    0.000    0.000    0.000 {method 'isidentifier' of 'str' objects}
			6    0.000    0.000    0.000    0.000 {method 'items' of 'dict' objects}
		  165    0.000    0.000    0.000    0.000 {method 'join' of 'str' objects}
		  405    0.000    0.000    0.000    0.000 {method 'lstrip' of 'str' objects}
		 2006    0.003    0.000    0.003    0.000 {method 'nonzero' of 'numpy.ndarray' objects}
		  116    0.000    0.000    0.000    0.000 {method 'pop' of 'dict' objects}
		  534    0.002    0.000    0.002    0.000 {method 'ravel' of 'numpy.ndarray' objects}
		 4800    4.358    0.001    4.358    0.001 {method 'reduce' of 'numpy.ufunc' objects}
		   18    0.000    0.000    0.000    0.000 {method 'replace' of 'str' objects}
		 2607    0.017    0.000    0.017    0.000 {method 'reshape' of 'numpy.ndarray' objects}
		 2000    0.006    0.000    0.006    0.000 {method 'sort' of 'numpy.ndarray' objects}
		  192    0.000    0.000    0.000    0.000 {method 'splitlines' of 'str' objects}
			6    0.000    0.000    0.000    0.000 {method 'take' of 'numpy.ndarray' objects}
		   45    0.000    0.000    0.000    0.000 {method 'values' of 'mappingproxy' objects}