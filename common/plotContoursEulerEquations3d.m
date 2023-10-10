function plotContoursEulerEquations3d(x,y,z,q,t)
    DIM = 3;
    colormap jet;
    subplot(231);
        isosurface(x,y,z,reshape(q(:,1),size(x)),0.3); hold on
        isosurface(x,y,z,reshape(q(:,1),size(x)),0.5);
        isosurface(x,y,z,reshape(q(:,1),size(x)),0.9); hold off
        view(DIM);
        title(sprintf('$\\rho(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        grid on
        clim([0.5,1.5]);
        colorbar;
    subplot(232);
        isosurface(x,y,z,reshape(q(:,2),size(x)),0.6);
        isosurface(x,y,z,reshape(q(:,2),size(x)),0.8);
        isosurface(x,y,z,reshape(q(:,2),size(x)),1.2);
        isosurface(x,y,z,reshape(q(:,2),size(x)),1.4);
        view(DIM);
        title(sprintf('$u(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        grid on
        clim([0,2]);
        colorbar;
    subplot(233);
        isosurface(x,y,z,reshape(q(:,3),size(x)),0.6);
        isosurface(x,y,z,reshape(q(:,3),size(x)),0.8);
        isosurface(x,y,z,reshape(q(:,3),size(x)),1.2);
        isosurface(x,y,z,reshape(q(:,3),size(x)),1.4);
        view(DIM);
        title(sprintf('$v(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        grid on
        clim([0,2]);
        colorbar;
    subplot(234);
        isosurface(x,y,z,reshape(q(:,4),size(x)),0.6);
        isosurface(x,y,z,reshape(q(:,4),size(x)),0.8);
        isosurface(x,y,z,reshape(q(:,4),size(x)),1.2);
        isosurface(x,y,z,reshape(q(:,4),size(x)),1.4);
        view(DIM);
        title(sprintf('$w(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        grid on
        clim([0,2]);
        colorbar;
    subplot(235);
        isosurface(x,y,z,reshape(q(:,5),size(x)),0.3);
        isosurface(x,y,z,reshape(q(:,5),size(x)),0.6);
        isosurface(x,y,z,reshape(q(:,5),size(x)),0.9);
        view(DIM);
        title(sprintf('$\\wp(x,y,t = %1.2f)$',t), Interpreter='latex');
        xlabel('$x$', Interpreter='latex'); 
        ylabel('$y$', Interpreter='latex');
        zlabel('$z$', Interpreter='latex');
        grid on
        clim([0,2]);
        colorbar;
end
