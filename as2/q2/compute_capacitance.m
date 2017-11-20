function C = compute_capacitance(file)
[Potential, S] = SIMPLE2D_M(file);

U = Potential(:, 4);
W = 0.5 * my_dot(my_transpose(U), my_dot(S, U));

e0 = 8.854e-12;
C = 2 * W * e0 / 15^2;

% Need to multiply by 4 since this is only a quarter of the system.
C = 4 * C;

return

function y = my_zeros(i_dim, j_dim)
  y = [];
  for i = 1:i_dim
    y = [y; 0];
  end

  if nargin == 2
    ys = [y];
    for j = 2:j_dim
      ys = [ys y];
    end
    y = ys;
  end
return

function y = my_transpose(x)
  [i_dim, j_dim] = size(x);
  y = my_zeros(j_dim, i_dim);
  for i = 1:i_dim
    for j = 1:j_dim
      y(j, i) = x(i, j);
    end
  end
return

function y = my_dot(a, b)
  [ai_dim, aj_dim] = size(a);
  [bi_dim, bj_dim] = size(b);
  y = my_zeros(ai_dim, bj_dim);
  for i = 1:ai_dim
    for j = 1:bj_dim
      x1 = a(i, :);
      x2 = my_transpose(b(:, j));
      y(i, j) = sum(x1 .* x2);
    end
  end
return