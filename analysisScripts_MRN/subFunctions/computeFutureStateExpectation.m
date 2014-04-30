function [expectationMatrix]=computeFutureStateExpectation(stateMatrix, pBlueBallOnNextDraw, pRedBallOnNextDraw)

%% expected state is a mixture of the + 1 red state and the + 1 blue state,
%% weighted by the respective probabilities of drawing these bead types.

expectationMatrix=stateMatrix(1:end-1, 2:end).*pRedBallOnNextDraw(1:end-1,1:end-1)+ ...
    stateMatrix(2:end,1:end-1).*pBlueBallOnNextDraw(1:end-1,1:end-1);
