function task1

    function kepler_elements(a, e, i, L, w, omega, mu, t)

        theta = omega * pi/180;
        g = (w - omega) * pi/180;
        i = i * pi/180;

        Theta = [ cos(theta) , -sin(theta) , 0   ;
                  sin(theta) ,  cos(theta) , 0   ;
                      0      ,      0      , 1  ];

        I = [ 1 ,   0    ,    0      ;
              0 , cos(i) , -sin(i)   ;
              0 , sin(i) ,  cos(i)  ]; 

        G = [ cos(g) , -sin(g) , 0   ;
              sin(g) ,  cos(g) , 0   ;
                0    ,    0    , 1  ];

        Q = Theta*I*G;
        
        gamma = 1 + mu;
        n = sqrt(gamma/a^3);
        to = ((w - L)/n)*pi/180;
        l = n*(t*2*pi - to);
        u = l + e*sin(l + e*sin(l + e*sin(l)));

        r = Q*a*[cos(u)-e ; sin(u)*sqrt(1-e^2) ; 0 ];
        v = Q*[-sin(u);cos(u)*sqrt(1-e^2);0]*a*n/(1-e*cos(u));

        disp('Coordinates (r)')
        disp(num2str(r))
        disp(['|r| = ', num2str(norm(r))])

        disp('Speed (v)')
        disp(num2str(v))
        disp(['|v| = ', num2str(norm(v))])
        
    end

    function calculations(data)

        time = 0.54757015742;
        
        kepler_elements(data(1), data(2), data(3), data(4), data(5), data(6), data(7), time)

    end

    d = [0.387  0.205 7.004  252.250 77.457  48.330  1/6023600;
         0.723  0.006 3.394  181.979 131.602 76.679  1/408523;
         1      0.016 0      100.464 102.937 0       1/328900.5;
         1.523  0.093 1.849  -4.553  -23.943 49.559  1/3098708;
         5.202  0.048 1.304  34.396  14.728  100.473 1/1047.34;
         9.536  0.053 2.485  49.954  92.598  113.662 1/3497.8;
         19.189 0.047 0.772  313.238 170.954 74.016  1/22902.9;
         30.069 0.008 1.770  -55.120 44.964  131.784 1/19402;
         39.482 0.248 17.140 238.929 224.068 110.303 1/135000000];

    planets = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'};

    for i = 1:9
        
        disp(char(planets(i))) 
        
        calculations(d(i,:))
        
    end

end