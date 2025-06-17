PortLineL = 3e-3;
PortLineW = 1.6e-3;

Whi = 0.2e-3;
Wlo = 3.2e-3;
Lstub = 6.35e-3;
Wstub = 0.238e-3;

L1  = 0.85e-3;
L2  = 3.22e-3;
L3  = 1.54e-3;
L4  = 3.39e-3;
L5  = 2.27e-3;
Lstub1 = 3.45e-3;


std_width = 10e-6;
std_thickness = 0e-6;
pstd_EpsilonR = 0.00;
n_parameters = 28;
freq = linspace(0.5e9,14e9,101);
for loop_idx = 1:1
    if loop_idx == 1
        n_samples = 1000;
        xi = lhsnorm(zeros(n_parameters,1), eye(n_parameters), n_samples);
    else
        n_samples = 1000;
        xi = randn(n_samples,n_parameters);
    end
    S11 = zeros(n_samples,length(freq));
    S21 = zeros(n_samples,length(freq));
    S12 = zeros(n_samples,length(freq));
    S22 = zeros(n_samples,length(freq));

    tic
    for sample_idx = 1:n_samples
        La = PortLineL + xi(sample_idx,1)*std_width;
        Wa = PortLineW + xi(sample_idx,2)*std_width;

        Lb = L1 + xi(sample_idx,3)*std_width;
        Wb = Whi + xi(sample_idx,4)*std_width;

        Lc = L2 + xi(sample_idx,5)*std_width;
        Wc = Wlo + xi(sample_idx,6)*std_width;

        Ll = Lstub1 + xi(sample_idx,7)*std_width;
        Wl = Wstub + xi(sample_idx,8)*std_width;

        Ld = L3 + xi(sample_idx,9)*std_width;
        Wd = Whi + xi(sample_idx,10)*std_width;

        Le = L4 + xi(sample_idx,11)*std_width;
        We = Wlo + xi(sample_idx,12)*std_width;

        Lf = L5 + xi(sample_idx,13)*std_width;
        Wf = Whi + xi(sample_idx,14)*std_width;

        Lg = L4 + xi(sample_idx,15)*std_width;
        Wg = Wlo + xi(sample_idx,16)*std_width;

        Lm = Lstub1 + xi(sample_idx,17)*std_width;
        Wm = Wstub + xi(sample_idx,18)*std_width;

        Lh = L3 + xi(sample_idx,19)*std_width;
        Wh = Whi + xi(sample_idx,20)*std_width;

        Li = L2 + xi(sample_idx,21)*std_width;
        Wi = Wlo + xi(sample_idx,22)*std_width;

        Lj = L1 + xi(sample_idx,23)*std_width;
        Wj = Whi + xi(sample_idx,24)*std_width;

        Lk = PortLineL + xi(sample_idx,25)*std_width;
        Wk = PortLineW + xi(sample_idx,26)*std_width;


        a  = traceRectangular('Length',La,'Width',Wa,'Center',[La/2,0]);
        b  = traceRectangular('Length',Lb,'Width',Wb,'Center',[La+Lb/2,0]);
        c  = traceRectangular('Length',Lc,'Width',Wlo,'Center',[La+Lb+Lc/2,0]);
        l = traceRectangular('Length',Wl,'Width',Ll,'Center',[La+Lb+Lc+Ld/2,Ll/2+Wd/2]);
        d  = traceRectangular('Length',Ld,'Width',Wd,'Center',[La+Lb+Lc+Ld/2,0]);
        e  = traceRectangular('Length',Le,'Width',We,'Center',[La+Lb+Lc+Ld+Le/2,0]);
        f  = traceRectangular('Length',Lf,'Width',Wf,'Center',[La+Lb+Lc+Ld+Le+Lf/2,0]);
        g  = traceRectangular('Length',Lg,'Width',Wg,'Center',[La+Lb+Lc+Ld+Le+Lf+Lg/2,0]);
        m = traceRectangular('Length',Wm,'Width',Lm,'Center',[La+Lb+Lc+Ld+Le+Lf+Lg+Lh/2,Lm/2+Wh/2]);
        h  = traceRectangular('Length',Lh,'Width',Wh,'Center',[La+Lb+Lc+Ld+Le+Lf+Lg+Lh/2,0]);
        i  = traceRectangular('Length',Li,'Width',Wi,'Center',[La+Lb+Lc+Ld+Le+Lf+Lg+Lh+Li/2,0]);
        j  = traceRectangular('Length',Lj,'Width',Wj,'Center',[La+Lb+Lc+Ld+Le+Lf+Lg+Lh+Li+Lj/2,0]);
        k  = traceRectangular('Length',Lk,'Width',Wk,'Center',[La+Lb+Lc+Ld+Le+Lf+Lg+Lh+Li+Lj+Lk/2,0]);
        compfiltShape2 = a+b+c+d+e+f+g+h+i+j+k+l+m;
        % figure;
        % show(compfiltShape2);

        compfilt2 = pcbComponent;
        d = dielectric('Teflon');
        d.EpsilonR  = 2.2*(1 + xi(sample_idx,27)*pstd_EpsilonR);
        d.Thickness = 0.508e-3 + xi(sample_idx,28)*std_thickness;


        GPL1 = La + Lb + Lc + Ld + Le + Lf + Lg + Lh + Li + Lj + Lk;
        GPW  = 20e-3;
        gnd  = traceRectangular('Length',GPL1,'Width',GPW/2,'Center',[GPL1/2,0]);
        compfilt2.BoardThickness = d.Thickness;
        compfilt2.Layers         = {compfiltShape2,d,gnd};
        compfilt2.BoardShape     = gnd;
        compfilt2.FeedDiameter   = Wk/2;
        compfilt2.FeedLocations  = [0,0,1,3;GPL1,0,1,3];
        compfilt2.ViaLocations   = [La+Lb+Lc+Ld/2,Ll-Wd/2,1,3; La+Lb+Lc+Ld+Le+Lf+Lg+Lh/2,Lm-Wh/2,1,3];
        compfilt2.ViaDiameter    = Wstub/2;
        % figure;
        % show(compfilt2);

        spar5 = sparameters(compfilt2,freq);
        % figure;
        % rfplot(spar5);
        S11(sample_idx,:) = reshape(spar5.Parameters(1,1,:),1,[]);
        S21(sample_idx,:) = reshape(spar5.Parameters(2,1,:),1,[]);
        S12(sample_idx,:) = reshape(spar5.Parameters(1,2,:),1,[]);
        S22(sample_idx,:) = reshape(spar5.Parameters(2,2,:),1,[]);
    end
    t1 = toc;
    sample_data.samples = xi;
    sample_data.S11 = S11;
    sample_data.S21 = S21;
    sample_data.S12 = S12;
    sample_data.S22 = S22;
    sample_data.sample_time = t1;

    if loop_idx == 1
        filename = 'Tests\Filter\data_training';
        save(filename,'-struct','sample_data');
    else
        filename = 'Tests\Filter\data_mc';
        save(filename,'-struct','sample_data');
    end
end