function [hr, hu, hv, hp] = plotEulerEquations2d(x,y,q,t)
    DIM = 2;
    colormap jet;
    subplot(221); 
        hr = surf(x,y,reshape(q(:,1),size(x)));
        view(DIM); 
        title(sprintf('$\\rho(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        shading interp;
        clim([0.5,1.5]);
        colorbar;
    subplot(222); 
        hu = surf(x,y,reshape(q(:,2),size(x)));
        view(DIM); 
        title(sprintf('$u(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        shading interp; 
        clim([0,2]);
        colorbar;
    subplot(223); 
        hv = surf(x,y,reshape(q(:,3),size(x)));
        view(DIM); 
        title(sprintf('$v(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        shading interp; 
        clim([0,2]);
        colorbar;
    subplot(224); 
        hp = surf(x,y,reshape(q(:,4),size(x)));
        view(DIM); 
        title(sprintf('$\\wp(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        shading interp; 
        clim([0,2]);
        colorbar;
end