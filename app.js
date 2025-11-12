/* app.js — JavaScript heat transfer calculator
   Runs fully in the browser, compatible with GitHub Pages.
*/

// Ensure DOM is loaded before attaching events
window.onload = function () {
  document.getElementById('runBtn').addEventListener('click', runCalc);
  document.getElementById('resetBtn').addEventListener('click', resetForm);
};

function resetForm() {
  document.getElementById('T_inf').value = 25;
  document.getElementById('T_des').value = 20;
  document.getElementById('Q_batt').value = 10;
  document.getElementById('num_batt').value = 100;
  document.getElementById('eps').value = 0.8;
  document.getElementById('Abs_paint').value = 0.6;
  document.getElementById('L').value = 6.0;
  document.getElementById('W').value = 2.44;
  document.getElementById('H').value = 2.59;
  
  document.getElementById('areasOut').textContent = '—';
  document.getElementById('tableWrap').innerHTML = '';
}

function runCalc() {
  // Input values
  const T_inf_C = parseFloat(document.getElementById('T_inf').value);
  const T_des_C = parseFloat(document.getElementById('T_des').value);
  const Q_batt = parseFloat(document.getElementById('Q_batt').value);
  const num_batt = parseFloat(document.getElementById('num_batt').value);
  const eps = parseFloat(document.getElementById('eps').value);
  const Abs_paint = parseFloat(document.getElementById('Abs_paint').value);
  const inch_to_m = 0.0254;
  const L = parseFloat(document.getElementById('L').value) * inch_to_m;
  const W = parseFloat(document.getElementById('W').value) * inch_to_m;
  const H = parseFloat(document.getElementById('H').value) * inch_to_m;

  const T_inf = T_inf_C + 273.15;
  const T_des = T_des_C + 273.15;
  const Q_tot_batt = Q_batt * num_batt;

  // Areas
  const A_wall1 = L * H;
  const A_wall2 = W * H;
  const A_roof = L * W;

  document.getElementById('areasOut').textContent =
    `A_wall1 = ${A_wall1.toFixed(3)} m²\n` +
    `A_wall2 = ${A_wall2.toFixed(3)} m²\n` +
    `A_roof  = ${A_roof.toFixed(3)} m²\n`;

  // Constants
  const R_poly = 1.2;
  const R_spr_foam = 1.2;
  const Cond_steel = 43;
  const Boltz = 5.67e-8;
  const h1 = 10;
  const h2 = 13;
  const delta = 0.002;
  const G = 1367;

  // Angle list
  const theta_deg = Array.from({ length: 19 }, (_, i) => i * 5);
  const theta_rad = theta_deg.map(d => d * Math.PI / 180);

  // Thermal resistance formulas
  const Rtot_no_ins = a => 1/(h1*a) + delta/(Cond_steel*a) + 1/(h2*a);
  const Rtot_ins_rf = a => 1/(h1*a) + delta/(Cond_steel*a) + R_spr_foam + 1/(h2*a);
  const Rtot_ins_wl = a => 1/(h1*a) + delta/(Cond_steel*a) + R_poly + 1/(h2*a);

  const Utot_no_ins_rf  = 1 / Rtot_no_ins(A_roof);
  const Utot_no_ins_wl1 = 1 / Rtot_no_ins(A_wall1);
  const Utot_no_ins_wl2 = 1 / Rtot_no_ins(A_wall2);
  const Utot_ins_rf  = 1 / Rtot_ins_rf(A_roof);
  const Utot_ins_wl1 = 1 / Rtot_ins_wl(A_wall1);
  const Utot_ins_wl2 = 1 / Rtot_ins_wl(A_wall2);

  // Simple bisection solver
  function solveBisection(f, a, b, tol=1e-6, maxIter=100) {
    let fa = f(a), fb = f(b);
    if (fa * fb > 0) return NaN;

    for (let i = 0; i < maxIter; i++) {
      const mid = (a + b) / 2;
      const fm = f(mid);
      if (Math.abs(fm) < tol) return mid;
      if (fa * fm <= 0) { b = mid; fb = fm; }
      else { a = mid; fa = fm; }
    }
    return (a + b) / 2;
  }

  // Solve surrounding temperature
  const T_surr_K1 = [];
  const T_surr_K2 = [];

  for (let i = 0; i < theta_rad.length; i++) {
    const θ = theta_rad[i];
    const Q_solar_roof = G * A_roof * Abs_paint * Math.sin(θ);
    const Q_solar_wall = G * A_wall1 * Abs_paint * Math.cos(θ);

    const f1 = Ts => eps * Boltz * A_roof * (Ts**4 - T_inf**4)
                      + h1*A_roof*(Ts - T_inf) - Q_solar_roof;

    const f2 = Ts => eps * Boltz * A_wall1 * (Ts**4 - T_inf**4)
                      + h1*A_wall1*(Ts - T_inf) - Q_solar_wall;

    T_surr_K1.push(solveBisection(f1, T_inf - 50, T_inf + 500));
    T_surr_K2.push(solveBisection(f2, T_inf - 50, T_inf + 500));
  }

  // Heat transfer arrays
  const Q_rf_no_ins = T_surr_K1.map(Ts => Utot_no_ins_rf * (Ts - T_des));
  const Q_rf_ins    = T_surr_K1.map(Ts => Utot_ins_rf    * (Ts - T_des));
  const Q_wl1_no_ins = T_surr_K2.map(Ts => Utot_no_ins_wl1 * (Ts - T_des));
  const Q_wl1_ins    = T_surr_K2.map(Ts => Utot_ins_wl1    * (Ts - T_des));

  const Q_wl2_no_ins = theta_deg.map(() => Utot_no_ins_wl2 * (T_inf - T_des));
  const Q_wl2_ins    = theta_deg.map(() => Utot_ins_wl2    * (T_inf - T_des));

  const Q_fl_no_ins = theta_deg.map(() => Utot_no_ins_rf * (T_inf - T_des));
  const Q_fl_ins    = theta_deg.map(() => Utot_ins_rf    * (T_inf - T_des));

  // Total Q
  const Q_tot_no_ins = theta_deg.map((_, i) =>
    Q_tot_batt + Q_rf_no_ins[i] + Q_wl1_no_ins[i] + Q_wl2_no_ins[i] + Q_fl_no_ins[i]
  );

  const Q_tot_ins = theta_deg.map((_, i) =>
    Q_tot_batt + Q_rf_ins[i] + Q_wl1_ins[i] + Q_wl2_ins[i] + Q_fl_ins[i]
  );

  // Build table
  let html = `<table>
    <thead>
      <tr>
        <th>Sun Angle (°)</th>
        <th>Q_no_ins (W)</th>
        <th>Q_ins (W)</th>
      </tr>
    </thead><tbody>`;

  for (let i = 0; i < theta_deg.length; i++) {
    html += `<tr>
      <td>${theta_deg[i]}</td>
      <td>${Q_tot_no_ins[i].toFixed(2)}</td>
      <td>${Q_tot_ins[i].toFixed(2)}</td>
    </tr>`;
  }

  html += "</tbody></table>";
  document.getElementById("tableWrap").innerHTML = html;
}


