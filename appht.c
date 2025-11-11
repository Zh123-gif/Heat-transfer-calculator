/* app.js
   Port of the MATLAB heat transfer calculation to JavaScript.
   Runs entirely in the browser. */
document.getElementById('runBtn').addEventListener('click', runCalc);
document.getElementById('resetBtn').addEventListener('click', resetForm);

function resetForm(){
  // reset to defaults used at load
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

function runCalc(){
  // Read inputs
  const T_inf_C = parseFloat(document.getElementById('T_inf').value);
  const T_des_C = parseFloat(document.getElementById('T_des').value);
  const Q_batt = parseFloat(document.getElementById('Q_batt').value);
  const num_batt = parseFloat(document.getElementById('num_batt').value);
  const eps = parseFloat(document.getElementById('eps').value);
  const Abs_paint = parseFloat(document.getElementById('Abs_paint').value);
  const L = parseFloat(document.getElementById('L').value);
  const W = parseFloat(document.getElementById('W').value);
  const H = parseFloat(document.getElementById('H').value);

  // Convert to Kelvin where needed (MATLAB added 273.15)
  const T_inf = T_inf_C + 273.15;
  const T_des = T_des_C + 273.15;

  const Q_tot_batt = Q_batt * num_batt;

  // Compute areas
  const A_wall1 = L * H;
  const A_wall2 = W * H;
  const A_roof  = L * W;

  document.getElementById('areasOut').textContent =
    `A_wall1 = ${A_wall1.toFixed(3)} m^2\n` +
    `A_wall2 = ${A_wall2.toFixed(3)} m^2\n` +
    `A_roof  = ${A_roof.toFixed(3)} m^2\n`;

  // Fixed parameters
  const R_poly = 1.2;
  const R_spr_foam = 1.2;
  const Cond_steel = 43;
  const Boltz = 5.67e-8;
  const h1 = 10;
  const h2 = 13;
  const delta = 0.002;
  const G = 1367;

  // Angles
  const theta_deg = Array.from({length: Math.floor(90/5)+1}, (_,i)=> i*5); // 0,5,...90
  const theta_rad = theta_deg.map(d => d * Math.PI/180);

  // Thermal resistances and U
  function Rtot_no_ins(area){
    return 1/(h1*area) + delta/Cond_steel/area + 1/(h2*area);
  }
  function Rtot_ins_rf(area){
    return 1/(h1*area) + delta/Cond_steel/area + R_spr_foam + 1/(h2*area);
  }
  function Rtot_ins_wl(area){
    return 1/(h1*area) + delta/Cond_steel/area + R_poly + 1/(h2*area);
  }

  const Rtot_no_ins_rf  = Rtot_no_ins(A_roof);
  const Rtot_no_ins_wl1 = Rtot_no_ins(A_wall1);
  const Rtot_no_ins_wl2 = Rtot_no_ins(A_wall2);

  const Rtot_ins_rf  = Rtot_ins_rf(A_roof);
  const Rtot_ins_wl1 = Rtot_ins_wl(A_wall1);
  const Rtot_ins_wl2 = Rtot_ins_wl(A_wall2);

  const Utot_no_ins_rf  = 1 / Rtot_no_ins_rf;
  const Utot_no_ins_wl1 = 1 / Rtot_no_ins_wl1;
  const Utot_no_ins_wl2 = 1 / Rtot_no_ins_wl2;

  const Utot_ins_rf  = 1 / Rtot_ins_rf;
  const Utot_ins_wl1 = 1 / Rtot_ins_wl1;
  const Utot_ins_wl2 = 1 / Rtot_ins_wl2;

  // Root solver (bisection) for energy balance: solves f(T_surr)=0
  function solveBisection(func, a, b, tol=1e-6, maxIter=200){
    let fa = func(a), fb = func(b);
    if (isNaN(fa) || isNaN(fb)) return NaN;
    // Ensure bracket; if not, expand a bit
    if (fa * fb > 0){
      // try expanding search window
      let factor = 1;
      for (let k=0;k<40;k++){
        factor *= 1.5;
        a = a - factor;
        b = b + factor;
        fa = func(a); fb = func(b);
        if (fa * fb <= 0) break;
      }
      if (fa * fb > 0) {
        // fallback: use secant starting near T_inf
        let x0 = a, x1 = b;
        for (let i=0;i<maxIter;i++){
          let f0 = func(x0), f1 = func(x1);
          if (Math.abs(f1 - f0) < 1e-12) break;
          let x2 = x1 - f1*(x1-x0)/(f1-f0);
          if (!isFinite(x2)) break;
          x0 = x1; x1 = x2;
          if (Math.abs(func(x1)) < tol) return x1;
        }
        return x1;
      }
    }
    let lo=a, hi=b;
    for (let i=0;i<maxIter;i++){
      let mid = 0.5*(lo+hi);
      let fm = func(mid);
      if (Math.abs(fm) < tol) return mid;
      if (fa * fm <= 0){
        hi = mid;
        fb = fm;
      } else {
        lo = mid;
        fa = fm;
      }
    }
    return 0.5*(lo+hi);
  }

  // Compute T_surr arrays
  const T_surr_K1 = new Array(theta_rad.length).fill(NaN);
  const T_surr_K2 = new Array(theta_rad.length).fill(NaN);

  for (let i=0;i<theta_rad.length;i++){
    const theta = theta_rad[i];

    const Q_solar_roof = G * A_roof  * Abs_paint * Math.sin(theta);
    const Q_solar_wall = G * A_wall1 * Abs_paint * Math.cos(theta);

    // Roof surrounding temperature energy balance
    const energy_balance1 = (T_surr) => eps*Boltz*A_roof*(Math.pow(T_surr,4) - Math.pow(T_inf,4)) + h1*A_roof*(T_surr - T_inf) - Q_solar_roof;
    // Wall surrounding temperature energy balance
    const energy_balance2 = (T_surr) => eps*Boltz*A_wall1*(Math.pow(T_surr,4) - Math.pow(T_inf,4)) + h1*A_wall1*(T_surr - T_inf) - Q_solar_wall;

    // initial bracket for bisection: [T_inf - 50, T_inf + 200] K
    const low = Math.max(1, T_inf - 200);
    const high = T_inf + 500;
    T_surr_K1[i] = solveBisection(energy_balance1, low, high);
    T_surr_K2[i] = solveBisection(energy_balance2, low, high);
  }

  // Heat transfers (vectorized analog)
  const Q_rf_no_ins = T_surr_K1.map(Ts => Utot_no_ins_rf  * (Ts - T_des));
  const Q_wl1a_no_ins = T_surr_K2.map(Ts => Utot_no_ins_wl1 * (Ts - T_des));
  const Q_wl2a_no_ins_ns = theta_deg.map(_=> Utot_no_ins_wl2 * (T_inf - T_des));
  const Q_wl1b_no_ins_ns = theta_deg.map(_=> Utot_no_ins_wl1 * (T_inf - T_des));
  const Q_wl2b_no_ins_ns = theta_deg.map(_=> Utot_no_ins_wl2 * (T_inf - T_des));
  const Q_fl_no_ins = theta_deg.map(_=> Utot_no_ins_rf * (T_inf - T_des));

  const Q_wl1a_ins = T_surr_K2.map(Ts => Utot_ins_wl1 * (Ts - T_des));
  const Q_wl2a_ins_ns = theta_deg.map(_=> Utot_ins_wl2 * (T_inf - T_des));
  const Q_wl1b_ins_ns = theta_deg.map(_=> Utot_ins_wl1 * (T_inf - T_des));
  const Q_wl2b_ins_ns = theta_deg.map(_=> Utot_ins_wl2 * (T_inf - T_des));

  const Q_rf_ins = T_surr_K1.map(Ts => Utot_ins_rf * (Ts - T_des));
  const Q_fl_ins = theta_deg.map(_=> Utot_ins_rf * (T_inf - T_des));

  // TOTALS
  const Q_tot_no_ins = theta_deg.map((_,i) =>
    Q_tot_batt + Q_wl2b_no_ins_ns[i] + Q_wl2a_no_ins_ns[i] +
    Q_rf_no_ins[i] + Q_wl1a_no_ins[i] + Q_wl1b_no_ins_ns[i] + Q_fl_no_ins[i]
  );

  const Q_tot_no_ins_rf = theta_deg.map((_,i) =>
    Q_tot_batt + Q_wl2b_ins_ns[i] + Q_wl2a_ins_ns[i] +
    Q_rf_no_ins[i] + Q_wl1a_ins[i] + Q_wl1b_ins_ns[i] + Q_fl_no_ins[i]
  );

  const Q_tot_ins = theta_deg.map((_,i) =>
    Q_tot_batt + Q_wl2b_ins_ns[i] + Q_wl2a_ins_ns[i] +
    Q_rf_ins[i] + Q_wl1a_ins[i] + Q_wl1b_ins_ns[i] + Q_fl_ins[i]
  );

  // Build HTML table
  let html = '<table><thead><tr>';
  html += '<th>Sun Angle (deg)</th><th>Q_no_ins (W)</th><th>Q_no_ins_roof (W)</th><th>Q_ins (W)</th>';
  html += '</tr></thead><tbody>';
  for (let i=0;i<theta_deg.length;i++){
    html += '<tr>';
    html += `<td style="text-align:center">${theta_deg[i]}</td>`;
    html += `<td>${Q_tot_no_ins[i].toFixed(2)}</td>`;
    html += `<td>${Q_tot_no_ins_rf[i].toFixed(2)}</td>`;
    html += `<td>${Q_tot_ins[i].toFixed(2)}</td>`;
    html += '</tr>';
  }
  html += '</tbody></table>';
  document.getElementById('tableWrap').innerHTML = html;
}