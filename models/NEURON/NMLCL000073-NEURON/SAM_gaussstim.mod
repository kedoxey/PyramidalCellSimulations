COMMENT
generate sinosoidal amplitude modulated pulse trains;

x=Xmax*[1+kcos(2*pi*fm(t-t_shift))], where x is time-varying envelope of a pulse train

Xmax=N(mean,std)

ENDCOMMENT


NEURON	{ 
  POINT_PROCESS SAM_GaussStim
  RANGE x
  RANGE interval, number, start,factor,refrac
  
  RANGE N_backward,N_forward,N_normal,N_total

  RANGE k,fm,t_shift
  RANGE amp_mean,amp_std, rand
  
}

PARAMETER {
	interval	= 10 (ms) <1e-9,1e9>: time between spikes (msec)
	number	= 10000 <0,1e9>	: number of spikes
	start		= 0 (ms)	: start of first spike
	factor =4  : portion of std to mean
	refrac		=0.5 <0,3>    : absolute refractory period, up limit of freq=1kHz
						: refrac >=dt, otherwise x can't be reset to 0 
        pi=3.1415927
	k=1		<0,1>   : modulation depth
	fm=200		<0, 1000>   : modulation frequency
	t_shift=0 (ms)
}

ASSIGNED {
	x
	event (ms)
	
	on
	end (ms)
	m (ms)            : mean of Gaussian
	diff (ms)
	N_forward  : swap spike whose value exceed forwardly one interval 
	N_backward  : swap spike whose value exceed backwardly one interval 
	N_normal
	N_total
	
 	amp_mean
	amp_std
	rand
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	x = 0
	diff=0
	m=interval/2  : each T has normal distribution N(interval/2, interval/2/factor)
	N_forward=0
	N_backward=0
	N_normal=0
	N_total=0
	if (start >= 0 && number > 0) {
		
		event = start + invl(m) 
		
		if (event < 0) {
			event = 0
		}
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = t
		end = t + 1e-6 + interval*(number-1)
		:printf("next event=%g\n",event)
	}
}

FUNCTION invl(mean (ms)) (ms) { LOCAL std
	if (mean <= 0.) {
		mean = .01 (ms) : I would worry if it were 0.
	}
	std=mean/factor  : std=T/(2*factor)
	invl = normrand(mean, std)  : relative to current interval 
	
	if(invl>=interval) { 
		:printf("original=%g\n",invl)
		invl=fmod(invl,interval)
		:printf("now=%g\n",invl)

		N_forward=N_forward+1


		}else if(invl<0) { 

			:printf("original=%g\n",invl)
			invl=fmod(invl,interval)+interval
			:printf("now=%g\n",invl)

			N_backward=N_backward+1
			}else {
			N_normal=N_normal+1
			}
		
		diff=interval-invl
	
	:printf("invl=%g\n",invl)
}

PROCEDURE event_time() {LOCAL diff2,T,rnd
        diff2=diff
	if (number > 0) {
	   T=invl(m)
	   rnd=T
	   if(T==0 && diff2==0) { T=T+dt } : previous and current spikes overlapped
	   
	   :printf("event=%g  diff=%g\n",event,diff2)
	   event = T+event + diff2    :compute absolute event time, relative to 0ms
	   
 	   N_total=N_total+1
 	}
 			
	if (event > end) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { : external event
		if (w > 0 && on == 0) { : turn on spike sequence
			init_sequence(t)
			net_send(0, 1)
		}else if (w < 0 && on == 1) { : turn off spiking
			on = 0
		}
	}
	if (flag == 3) { : from INITIAL
		if (on == 0) {
			init_sequence(t)
			net_send(0, 1)
		}
	}
	if (flag == 1 && on == 1) {
		if(x == 0){  : after refractory
				rand=normrand(amp_mean, amp_std) : normal distribution
				x=rand*(1+k*cos((0.001)*2*pi*fm*(t-t_shift)))/(1+k)
				net_event(t)
				event_time()  : after each spike call next spike time

		  		 if (event-t <= refrac+dt && event-t >= refrac) {
		      		 	event=event+dt    : happen when next event at reset edge time
		       			:printf("next spike on edge %g\n", event)
		  		 }

		  		 if (on==1) {
						net_send(event - t, 1) : self-feed next event
		  				 }
				net_send(refrac, 2)  : refractory period
		
		} else if (x!=0) {
			net_event(t)
			event_time() : although this spike will be ignored, still call next spike : independent cycles 
			:printf("inside refrac\n")

			if (on == 1) {
				net_send(event - t, 1)
					  }
		}
		:printf("x=1 at %g\n",t)
	}
	if (flag == 2) {
		:printf("x=0 at %g\n",t)
		x = 0
		
	}
}



