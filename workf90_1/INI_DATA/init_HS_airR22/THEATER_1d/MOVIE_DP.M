% program movie_dp

% Record the movie
for jjj = 1:5
    show_dp;
    MM(jjj) = getframe;
    pause;
end
% Play the movie twenty times
movie(MM,3,4);

movie2avi(MM,'mvp102','fps',4)