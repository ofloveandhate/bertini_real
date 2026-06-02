i have played with the decomposition a lot.  it's broken.  the first noted failure on verb level 0 is from curve interslicing on the sphere

```
connecting midpoint downstairs, 0 of 13
num_edges = 35
done decomposing critical curve
interslicing sphere
curve_interslice_crit_downstairs = [...
		-2.2914391e3+1i*0;
		-2.0264151e3+1i*0;
		-1.3853139e-1+1i*0;
		0+1i*0;
		3.8836141e-1+1i*0;
		6.2142549e1+1i*0;
		6.2149054e1+1i*0;
		7.2639376e1+1i*0;
		4.3606150e2+1i*0;
		4.3659728e2+1i*0;
		2.0629506e3+1i*0;
		2.0647338e3+1i*0;
		2.0697654e3+1i*0;
		2.2913987e3+1i*0;
];

there were non-unique midslice points in interval 5.
trying to recover the failure by tightening tracking tolerances...
new temporary tracktolBEFOREeg: 1e-12 tracktolDURINGeg: 1e-14
there were non-unique midslice points in interval 5.  your decomposition is possibly incorrect about the missed points, if the path crossings obscured real points
connecting midpoint downstairs, 0 of 13
num_edges = 94
```




here's the verb level 4 output for the inter

```
solution 0, success 1, multi 2, isFinite 1, isSing 1, isReal 0, cycle_num 2
solution 1, success 1, multi 1, isFinite 1, isSing 0, isReal 1, cycle_num 1
solution 2, success 1, multi -1, isFinite 1, isSing 1, isReal 0, cycle_num 2
solution 3, success 1, multi 1, isFinite 1, isSing 0, isReal 1, cycle_num 1
solution 4, success 1, multi 1, isFinite 1, isSing 0, isReal 0, cycle_num 1
solution 5, success 1, multi 1, isFinite 1, isSing 0, isReal 0, cycle_num 1
witness set has 4 total variables, 4 natural variables.
******
6 points
******
point_0 = [...
		3.5905892e3+1i*1.6107098e-6;
		2.8614347e2+1i*-5.3138503e-5;
		1.5746248e2+1i*5.8624102e-5;
];

point_1 = [...
		3.5385138e3+1i*0;
		-1.4180086e2+1i*0;
		6.2270744e2+1i*-3.9806008e-37;
];

point_2 = [...
		3.5905892e3+1i*1.6107098e-6;
		2.8614347e2+1i*-5.3138503e-5;
		1.5746248e2+1i*5.8624102e-5;
];


```




getting this:

```

*****************************
midslice 5 / 13, edge 0 / 18
current midpoint: 177 
bottom system is "input_midslice_5", which is not in md_config
	top system: input_surf_sphere
	bottom system: input_midslice_5
0 bottom variables
bailing out 5 0.
tracking from these point indices:
174 177 55
```

midslice 5 edge 0 is bad.  top is surface, bottom is just itself again.  whis is point 174.  174.  




getting closer.  i found that a point should have been added elsewhere


Out[14]: 
['input_critical_curve',
 'input_singcurve_mult_2_0',
 'input_surf_sphere',
 'W_total_crit_nonexistant_filename',
 'should_have_already_been_added_elsewhere',
 'input_midslice_0',
 'input_midslice_1',
 'input_midslice_2',
 'input_midslice_3',
 'input_midslice_4',
 'input_midslice_5',




