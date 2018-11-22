
function qdot = qdotfunc(q, u, r)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    s1 = sin(q(1));
    s2 = sin(q(2));
    s3 = sin(q(3));

    c1 = cos(q(1));
    c2 = cos(q(2));
    c3 = cos(q(3));

    A = [c3/c2, -s3/c2, 0, 0, 0;
         s3,     c3,    0, 0, 0;
         -c3*s2/c2, s3*s2/c2, 1, 0, 0;
         r*(s1*s3-c1*s2*c3), r*(s1*c3+c1*s2*s3), r*c1*c2, 1, 0;
         -r*(c1*s3+s1*s2*c3), -r*(c1*c3-s1*s2*s3), r*(s1*c2), 0, 1];
     
     qdot = A*u;

end