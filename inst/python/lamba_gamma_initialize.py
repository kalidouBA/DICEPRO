from hyperopt import hp
from hyperopt.pyll.stochastic import sample
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio

def custom_space():
    gamma = hp.loguniform("gamma", np.log(1), np.log(1e5))
    lambda_factor = hp.loguniform("lambda_factor", np.log(2), np.log(1e2))  # lambda_ = gamma * factor
    lambda_ = gamma * lambda_factor
    p_prime = hp.loguniform("p_prime", np.log(1e-1), np.log(1))
    return {
        "lambda_": lambda_,
        "gamma": gamma,
        "p_prime": p_prime
      }  
    
    
space = custom_space()
samples = [sample(space) for _ in range(200)]
gamma_samples = [s["gamma"] for s in samples]
lambda_samples = [s["lambda_"] for s in samples]


gamma_range = np.logspace(0, 5, 100)
lambda_min = 2 * gamma_range
lambda_max = 1e2 * gamma_range


scatter = go.Scatter(
    x=gamma_samples,
    y=lambda_samples,
    mode='markers',
    marker=dict(size=3, color='royalblue', opacity=0.4),
    name='Samples'
)

line_min = go.Scatter(
    x=gamma_range,
    y=lambda_min,
    mode='lines',
    line=dict(dash='dash', color='red'),
    name='λ = 2 × γ'
)

line_max = go.Scatter(
    x=gamma_range,
    y=lambda_max,
    mode='lines',
    line=dict(dash='dash', color='green'),
    name='λ = 10⁴ × γ'
)

layout = go.Layout(
    title="Échantillons (γ, λ) avec contrainte λ = k·γ, k ∈ [10, 10⁴]",
    xaxis=dict(title='γ', type='log'),
    yaxis=dict(title='λ', type='log'),
    hovermode='closest'
)

fig = go.Figure(data=[scatter, line_min, line_max], layout=layout)

pio.write_html(fig, file="gamma_lambda_interactive_plot.html", auto_open=True)
