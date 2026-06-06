// Alkaline Water Electrolysis — Electrochemical Model
// Based on published literature (LeRoy et al. 1980, Zeng & Zhang 2010)
// Membrane: Zirfon PERL UTP 500 (~0.46 mm); Electrolyte: 30 wt% KOH

const Rc=8.314, Fc=96485, Tref=298.15, Pref=1.01325, Eth=1.48, Ae=0.0729, Nc=20, dm=0.46e-3;

function calc(j, Tc, Pb) {
  const T = Tc + 273.15;
  const m = (30*(183.1221 - 0.56845*T + 984.5679*Math.exp(30/115.96277))) / 5610.5;
  const bub = Math.max(0, 0.023 * Math.pow(Math.max(j*1e4,0.01),0.3)
              * Math.pow((T/Tref)*(Pref/Pb), 2/3));
  const sKOH = (-2.041*m - 2.8e-3*m*m + 5.332e-3*m*T + 207.2*m/T
               + 1.043e-3*m*m*m - 3e-7*m*m*T*T) * 100;
  const sE = 6e7 - 279650*T + 532*T*T - 0.38057*T*T*T;
  const joc = 1.5e-4 * Math.pow(Pb/Pref,0.1) * Math.exp(-23000/(Rc*T)*(1-T/Tref));
  const joa = 9e-5  * Math.pow(Pb/Pref,0.1) * Math.exp(-42000/(Rc*T)*(1-T/Tref));
  const eta_c = (Rc*T/(0.5*Fc)) * Math.log(Math.max(j/(1-bub), joc*1.001) / joc);
  const eta_a = (Rc*T/(0.5*Fc)) * Math.log(Math.max(j/(1-bub), joa*1.001) / joa);
  const act_ov = eta_c + eta_a;
  const R_KOH = dm / (sKOH * Math.pow(Math.max(1-bub,0.01), 1.5));
  const R_elec = dm / Math.max(sE, 1);
  const ohm_ov = (j*1e4) * (R_KOH + R_elec);
  const PH2O = Math.pow(10, 8.07131 - 1730.63/(233.426+Tc)) * 0.00133322;
  const Er0  = 1.50342 - 9.956e-4*T + 2.5e-7*T*T;
  const Erev = Er0 + (Rc*T/(2*Fc)) *
    Math.log(Math.max(Math.pow(Pb-PH2O,2)*Math.pow(Pb-PH2O,0.5)/Math.max(PH2O,1e-10), 1e-10));
  const V_cell = Erev + act_ov + ohm_ov;
  const I = (j*1e4)*Ae, V_stack = Nc*V_cell, P_stack_kW = V_stack*I/1000;
  const h2_mol=I/(2*Fc), h2_kg_hr=h2_mol*0.002016*3600, h2_Nm3=h2_mol*0.022414*3600;
  const SEC = P_stack_kW / h2_kg_hr;
  const eff = (Eth / V_cell) * 100;
  return {bub, m, sKOH, sE, eta_c, eta_a, act_ov, ohm_ov, Erev, V_cell,
          I, V_stack, P_stack_kW, h2_mol, h2_kg_hr, h2_Nm3, SEC, eff};
}

let T=80, P=3, j=0.10, isPaused=false;
const HIST={t:[],V:[],E:[],H:[],Pk:[],S:[],B:[]}, MAXH=60;
const liveCh={}, swCh={}, sensCh={};
let fcCh=null;

function noise(s) { return (Math.random()-0.5)*2*s; }

function togglePause() {
  isPaused = !isPaused;
  const btn=document.getElementById('pause_btn'),
        dot=document.getElementById('pulse_dot'),
        lbl=document.getElementById('live_lbl'),
        ban=document.getElementById('paused_banner');
  if (isPaused) {
    btn.textContent='▶ Resume'; btn.classList.add('paused');
    dot.style.animation='none'; dot.style.background='#f59e0b';
    lbl.textContent='PAUSED'; ban.classList.add('visible');
  } else {
    btn.textContent='⏸ Pause'; btn.classList.remove('paused');
    dot.style.animation='pulse 1.5s infinite'; dot.style.background='#10b981';
    lbl.textContent='LIVE'; ban.classList.remove('visible');
    updateUI();
  }
}

function mkLiveChart(id, lbl, col) {
  const ctx = document.getElementById(id); if(!ctx) return null;
  return new Chart(ctx.getContext('2d'), {
    type:'line',
    data:{labels:[], datasets:[{label:lbl, data:[], borderColor:col,
      backgroundColor:col+'20', borderWidth:1.5, pointRadius:0, fill:true, tension:0.4}]},
    options:{animation:{duration:0}, responsive:true, maintainAspectRatio:true,
      plugins:{legend:{display:false}},
      scales:{x:{display:false},
        y:{grid:{color:'rgba(255,255,255,0.04)'},
          ticks:{font:{size:9}, color:'#64748b', maxTicksLimit:4}}}}
  });
}
liveCh.V  = mkLiveChart('ch_V',  'V_cell','#3b82f6');
liveCh.E  = mkLiveChart('ch_E',  'Eff',   '#10b981');
liveCh.H  = mkLiveChart('ch_H',  'H₂',    '#06b6d4');
liveCh.Pk = mkLiveChart('ch_P2', 'Power', '#8b5cf6');
liveCh.S  = mkLiveChart('ch_S',  'SEC',   '#f59e0b');
liveCh.B  = mkLiveChart('ch_B',  'Bub',   '#ef4444');

function st(id,v){const e=document.getElementById(id);if(e)e.textContent=v;}
function sh(id,v){const e=document.getElementById(id);if(e)e.innerHTML=v;}
function badge(id,cls,txt){sh(id,`<span class="kpi-badge ${cls}">${txt}</span>`);}

function onSlide() {
  T=+document.getElementById('sl_T').value;
  P=+document.getElementById('sl_P').value;
  j=+document.getElementById('sl_j').value;
  st('T_disp', T+' °C');
  st('P_disp', P.toFixed(1)+' bar');
  st('j_disp', j.toFixed(2)+' A/cm²');
  updateUI();
}

function updateUI() {
  if (isPaused) return;
  const r = calc(j,T,P);
  const ts = new Date().toLocaleTimeString();
  st('clock', ts);

  const h2d = Math.max(0, r.h2_kg_hr + noise(r.h2_kg_hr*0.0018));
  const efd = Math.max(0, r.eff + noise(0.04));

  sh('k_Vcell', r.V_cell.toFixed(4)+'<span class="kpi-unit">V</span>');
  sh('k_eff',   efd.toFixed(2)+'<span class="kpi-unit">%</span>');
  sh('k_h2',    h2d.toFixed(5)+'<span class="kpi-unit">kg/hr</span>');
  sh('k_P',     r.P_stack_kW.toFixed(2)+'<span class="kpi-unit">kW</span>');
  sh('k_SEC',   r.SEC.toFixed(0)+'<span class="kpi-unit">kWh/kg</span>');
  sh('k_act',   r.act_ov.toFixed(4)+'<span class="kpi-unit">V</span>');
  sh('k_ohm',   r.ohm_ov.toFixed(5)+'<span class="kpi-unit">V</span>');
  sh('k_hmol',  (r.h2_mol*1000).toFixed(4)+'<span class="kpi-unit">mmol/s</span>');

  badge('kb_eff', efd>=80?'bg-good':efd>=70?'bg-warn':'bg-bad',
    efd>=80?'✓ Excellent':efd>=70?'⚠ OK':'↓ Low');
  badge('kb_V', r.V_cell<=1.85?'bg-good':r.V_cell<=2.1?'bg-warn':'bg-bad',
    r.V_cell<=1.85?'✓ Low':r.V_cell<=2.1?'⚠ Med':'↑ High');
  badge('kb_h2','bg-good', h2d>0.005?'↑ High':'→ Norm');
  badge('kb_P', 'bg-warn', r.P_stack_kW.toFixed(1)+' kW');
  badge('kb_SEC', r.SEC<950?'bg-good':r.SEC<1100?'bg-warn':'bg-bad',
    r.SEC<950?'✓ Low':r.SEC<1100?'⚠ Med':'↑ High');

  st('v_Erev',  r.Erev.toFixed(4)+' V');
  st('v_etac',  r.eta_c.toFixed(4)+' V');
  st('v_etaa',  r.eta_a.toFixed(4)+' V');
  st('v_ohm',   r.ohm_ov.toFixed(5)+' V');
  st('v_Vcell', r.V_cell.toFixed(4)+' V');

  const pE=(r.Erev/r.V_cell*100).toFixed(1),
        pA=(r.act_ov/r.V_cell*100).toFixed(1),
        pO=(r.ohm_ov/r.V_cell*100).toFixed(1);
  sh('vbar_el',
    `<div class="vbar-s" style="width:${pE}%;background:#10b981">${pE}%</div>
     <div class="vbar-s" style="width:${pA}%;background:#ef4444">${pA}%</div>
     <div class="vbar-s" style="width:${pO}%;background:#3b82f6">${pO}%</div>`);

  st('s_bub',   (r.bub*100).toFixed(3)+'%');
  st('s_sKOH',  r.sKOH.toFixed(1)+' S/m');
  st('s_Icell', r.I.toFixed(2)+' A');
  st('s_Vstack',r.V_stack.toFixed(3)+' V');
  st('s_Nm3',   r.h2_Nm3.toFixed(5)+' Nm³/hr');
  st('s_hmol',  (r.h2_mol*1e6).toFixed(3)+' μmol/s');

  const rf = calc(0.10,80,3);
  const pp = [['V_cell',rf.V_cell.toFixed(4),'1.732'],['Eff%',rf.eff.toFixed(2),'81.34'],
              ['H₂',rf.h2_kg_hr.toFixed(5),'0.00274'],['Power',rf.P_stack_kW.toFixed(2),'10.10'],
              ['E_rev',rf.Erev.toFixed(4),'~1.22']];
  sh('cmp_body', pp.map(([p,mv,pv])=>{
    const np=parseFloat(pv),nm=parseFloat(mv);
    const err=isNaN(np)?'—':(Math.abs(nm-np)/Math.abs(np)*100).toFixed(1)+'%';
    const ne=parseFloat(err);
    const cl=isNaN(ne)?'':ne<=2?'eg':ne<=5?'ew':'eb';
    return `<tr><td>${p}</td><td style="color:var(--txt1);text-align:center">${mv}</td>
            <td style="text-align:center">${pv}</td><td class="${cl}" style="text-align:right">${err}</td></tr>`;
  }).join(''));

  HIST.t.push(ts); HIST.V.push(r.V_cell); HIST.E.push(efd);
  HIST.H.push(h2d); HIST.Pk.push(r.P_stack_kW); HIST.S.push(r.SEC); HIST.B.push(r.bub);
  if(HIST.t.length>MAXH) Object.keys(HIST).forEach(k=>HIST[k].shift());
  ['V','E','H','Pk','S','B'].forEach(k=>{
    if(liveCh[k]){liveCh[k].data.labels=HIST.t;liveCh[k].data.datasets[0].data=HIST[k];liveCh[k].update('none');}
  });

  mlPredict(j,T,P,r);
  checkAnom(r);
  updateForecast(r);
  updateSens(j,T,P);
  st('st_upd','Last update: '+ts);
}

function showTab(name,el){
  document.querySelectorAll('.tab').forEach(t=>t.classList.remove('active'));
  document.querySelectorAll('.tab-pane').forEach(p=>p.classList.remove('active'));
  el.classList.add('active');
  document.getElementById('tab_'+name).classList.add('active');
  if(name==='sweep') setTimeout(runSweep,80);
  if(name==='sens')  setTimeout(()=>updateSens(j,T,P),80);
}

const ylabels={V_cell:'Cell Voltage (V)',eff:'Efficiency (%)',h2_kg_hr:'H₂ (kg/hr)',
  P_stack_kW:'Stack Power (kW)',SEC:'SEC (kWh/kg)',act_ov:'Act.OV (V)',ohm_ov:'Ohmic OV (V)'};

function mkLC2(id,dsets,xl,yl){
  const c=document.getElementById(id); if(!c)return;
  if(swCh[id]) swCh[id].destroy();
  swCh[id]=new Chart(c.getContext('2d'),{type:'line',
    data:{labels:dsets[0].labels, datasets:dsets.map(d=>({label:d.label,data:d.data,
      borderColor:d.c,backgroundColor:'transparent',borderWidth:1.5,pointRadius:0,tension:0.3}))},
    options:{responsive:true,maintainAspectRatio:true,
      plugins:{legend:{display:true,labels:{font:{size:9},color:'#94a3b8',boxWidth:10}}},
      scales:{x:{title:{display:true,text:xl,font:{size:10},color:'#64748b'},
        ticks:{maxTicksLimit:8,font:{size:8},color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}},
              y:{title:{display:true,text:yl,font:{size:10},color:'#64748b'},
        ticks:{font:{size:8},color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}}}}});
}

function runSweep(){
  const xp=document.getElementById('sw_x').value, yp=document.getElementById('sw_y').value;
  let xs=[],ys=[],xl='';
  if(xp==='j'){for(let jj=0.05;jj<=1.0;jj+=0.02){xs.push(+jj.toFixed(2));ys.push(calc(jj,T,P)[yp]);}xl='j (A/cm²)';}
  else if(xp==='T'){for(let tt=20;tt<=90;tt++){xs.push(tt);ys.push(calc(j,tt,P)[yp]);}xl='T (°C)';}
  else{for(let pp=1;pp<=10;pp+=0.2){xs.push(+pp.toFixed(1));ys.push(calc(j,T,pp)[yp]);}xl='P (bar)';}
  const c=document.getElementById('ch_sweep'); if(!c)return;
  if(swCh['main']) swCh['main'].destroy();
  swCh['main']=new Chart(c.getContext('2d'),{type:'line',
    data:{labels:xs, datasets:[{label:ylabels[yp],data:ys,borderColor:'#3b82f6',
      backgroundColor:'#3b82f618',borderWidth:2,pointRadius:0,fill:true,tension:0.3}]},
    options:{responsive:true,maintainAspectRatio:true,
      plugins:{legend:{display:true,labels:{color:'#94a3b8'}}},
      scales:{x:{title:{display:true,text:xl,color:'#64748b'},ticks:{color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}},
              y:{title:{display:true,text:ylabels[yp],color:'#64748b'},ticks:{color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}}}}});
  const jA=[]; for(let jj=0.05;jj<=1.0;jj+=0.02) jA.push(+jj.toFixed(2));
  const Ts=[40,60,80], cs=['#06b6d4','#f59e0b','#10b981'];
  mkLC2('ch_sw1',Ts.map((t,i)=>({label:`T=${t}°C`,data:jA.map(jj=>calc(jj,t,P).V_cell),c:cs[i],labels:jA})),'j (A/cm²)','V_cell (V)');
  mkLC2('ch_sw2',Ts.map((t,i)=>({label:`T=${t}°C`,data:jA.map(jj=>calc(jj,t,P).eff),  c:cs[i],labels:jA})),'j (A/cm²)','Eff (%)');
  const Ps=[1,3,6];
  mkLC2('ch_sw3',Ps.map((p,i)=>({label:`P=${p}bar`,data:jA.map(jj=>calc(jj,T,p).h2_kg_hr),c:cs[i],labels:jA})),'j (A/cm²)','H₂ (kg/hr)');
}

function mlPredict(j,T,P,r){
  const Vp = 1.1896 + 0.3124*j + 0.1872*j*j - 0.001483*(T-60)
            - 0.0098*Math.log(P+1) + 0.00021*(T-60)*j;
  const ep  = (Eth/Math.max(Vp,1.2))*100;
  const h2p = r.h2_mol*0.002016*3600;
  const errV = (Math.abs(Vp-r.V_cell)/r.V_cell*100);
  st('ml_V',   Vp.toFixed(4)+' V  (err: '+errV.toFixed(2)+'%)');
  st('ml_eff', ep.toFixed(2)+'%');
  st('ml_h2',  h2p.toFixed(5)+' kg/hr');
}

function runOpt(){
  let best={eff:0};
  for(let t=55;t<=85;t+=5)
    for(let p=2;p<=8;p++)
      for(let jj=0.05;jj<=0.30;jj+=0.025){
        const r=calc(jj,t,p);
        if(r.eff>best.eff && r.V_cell<2.1)
          best={eff:r.eff,j:jj,T:t,P:p,V:r.V_cell,SEC:r.SEC};
      }
  sh('opt_out',
    `<b style="color:#34d399">Optimal Found:</b><br>
    T = <b style="color:#f59e0b">${best.T}°C</b> · P = <b style="color:#8b5cf6">${best.P}bar</b>
    · j = <b style="color:#60a5fa">${best.j.toFixed(3)}A/cm²</b><br>
    η = <b style="color:#34d399">${best.eff.toFixed(2)}%</b>
    · V = <b>${best.V.toFixed(4)}V</b><br>
    SEC = <b style="color:#f59e0b">${best.SEC.toFixed(0)} kWh/kg</b>`);
}

function checkAnom(r){
  const a=[];
  if(r.eff<65)     a.push('<span style="color:#f87171">⚠ Efficiency below 65%</span>');
  if(r.V_cell>2.3) a.push('<span style="color:#f87171">⚠ Cell voltage >2.3V</span>');
  if(r.bub>0.4)    a.push('<span style="color:#fbbf24">⚠ Bubble coverage >40%</span>');
  if(r.act_ov>1.0) a.push('<span style="color:#fbbf24">⚠ Act. OV >1.0V</span>');
  if(r.eff>80)     a.push('<span style="color:#34d399">✅ Efficiency optimal</span>');
  if(r.V_cell<=1.85) a.push('<span style="color:#34d399">✅ Cell voltage nominal</span>');
  if(a.length===0) a.push('<span style="color:#34d399">✅ All parameters nominal</span>');
  sh('anom_out', a.join('<br>'));
}

async function evaluatePerformance(q){
  const btn=document.getElementById('report_btn'), out=document.getElementById('report_out');
  const r=calc(j,T,P);
  btn.disabled=true; btn.textContent='Generating...';
  out.textContent='Generating performance report...';
  const analysisConfig=`Electrolyzer performance evaluation module.
η=E_th/V_cell×100, E_th=1.48V. Base paper: j=0.1A/cm², T=80°C → η≈82.2%, V≈1.75V.
Generate performance observations and recommendations.`;
  const analysisInput= q
    ? `T=${T}°C P=${P}bar j=${j}A/cm² | V=${r.V_cell.toFixed(4)}V η=${r.eff.toFixed(2)}%
       SEC=${r.SEC.toFixed(0)} H₂=${r.h2_kg_hr.toFixed(5)}kg/hr
       η_act=${r.act_ov.toFixed(4)}V η_ohm=${r.ohm_ov.toFixed(5)}V θ=${(r.bub*100).toFixed(2)}%
       Q: ${q}`
    : `Analyze: T=${T}°C P=${P}bar j=${j}A/cm² V=${r.V_cell.toFixed(4)}V η=${r.eff.toFixed(2)}%
       SEC=${r.SEC.toFixed(0)} H₂=${r.h2_kg_hr.toFixed(5)}kg/hr
       η_act=${r.act_ov.toFixed(4)}V(${(r.act_ov/r.V_cell*100).toFixed(1)}%)
       η_ohm=${r.ohm_ov.toFixed(5)}V θ=${(r.bub*100).toFixed(2)}%.
       Dominant losses, performance quality, recommendations.`;
  try {
    const response=await fetch('/api/performance-report',{method:'POST',
      headers:{'Content-Type':'application/json'},
      body:JSON.stringify({analysisType:'performance-report',reportLength:600,
        configuration:analysisConfig,operatingConditions:analysisInput})});
    const result=await response.json();
    out.textContent=result.content?.[0]?.text||'Analysis unavailable.';
  } catch(e) {
    out.textContent='[Performance Evaluation]\n'
      +'• Efficiency: '+r.eff.toFixed(2)+'% — '+(r.eff>80?'Excellent, matches base paper ~82.2%':'Below target')+'\n'
      +'• Dominant loss: '+(r.act_ov>r.ohm_ov?'Activation OV '+r.act_ov.toFixed(4)+'V = '+(r.act_ov/r.V_cell*100).toFixed(1)+'%':'Ohmic OV')+'\n'
      +'• Bubble coverage: '+(r.bub*100).toFixed(2)+'% '+(r.bub>0.2?'(significant blockage)':'(acceptable)')+'\n'
      +'• Recommendation: '+(j<0.15?'Low j → high η. Increase j for higher H₂ output.':'Moderate j. Raise T to improve η.');
  }
  btn.disabled=false; btn.textContent='▶ Generate Report';
  document.getElementById('report_query').value='';
}

function updateSens(j0,T0,P0){
  const r0=calc(j0,T0,P0), h0=r0.h2_kg_hr;
  const dj=j0*0.01, dT=0.5, dP=P0*0.01;
  const dhj=(calc(j0+dj,T0,P0).h2_kg_hr-calc(j0-dj,T0,P0).h2_kg_hr)/(2*dj);
  const dhT=(calc(j0,T0+dT,P0).h2_kg_hr-calc(j0,T0-dT,P0).h2_kg_hr)/(2*dT);
  const dhP=(calc(j0,T0,P0+dP).h2_kg_hr-calc(j0,T0,P0-dP).h2_kg_hr)/(2*dP);
  const Sj=Math.abs(dhj)*j0/h0, ST=Math.abs(dhT)*T0/h0, SP=Math.abs(dhP)*P0/h0;
  const tot=Sj+ST+SP||1;
  const rows=[
    {l:'Current Density j',v:Sj,p:Sj/tot*100,d:dhj,c:'#3b82f6'},
    {l:'Temperature T',    v:ST,p:ST/tot*100,d:dhT,c:'#f59e0b'},
    {l:'Pressure P',       v:SP,p:SP/tot*100,d:dhP,c:'#8b5cf6'}
  ].sort((a,b)=>b.v-a.v);

  sh('sens_bars', rows.map(r=>
    `<div class="sens-row">
      <div class="sens-lbl">${r.l}</div>
      <div class="sens-track"><div class="sens-fill" style="width:${r.p.toFixed(0)}%;background:${r.c}"></div></div>
      <div class="sens-pct" style="color:${r.c}">${r.p.toFixed(1)}%</div>
      <div class="sens-val">S=${r.v.toFixed(3)}</div>
    </div>`).join(''));

  sh('sobol_table',
    `<table class="vtbl"><tr><th>Parameter</th><th style="text-align:center">Elasticity</th>
     <th style="text-align:right">∂ṁH₂/∂xᵢ</th><th style="text-align:center">Rank</th></tr>
     ${rows.map((r,i)=>`<tr>
       <td style="color:${r.c};font-weight:700">${r.l}</td>
       <td style="text-align:center;color:#e2e8f0">${r.v.toFixed(4)}</td>
       <td style="text-align:right;font-family:monospace;color:#94a3b8">${r.d.toExponential(3)}</td>
       <td style="text-align:center">${i+1}</td></tr>`).join('')}</table>`);

  function mkSC(id,xs,ys,xl,yl,col){
    const c=document.getElementById(id); if(!c)return;
    if(sensCh[id]) sensCh[id].destroy();
    sensCh[id]=new Chart(c.getContext('2d'),{type:'line',
      data:{labels:xs, datasets:[{data:ys,borderColor:col,backgroundColor:col+'18',
        borderWidth:1.5,pointRadius:0,fill:true,tension:0.3}]},
      options:{responsive:true,maintainAspectRatio:true,plugins:{legend:{display:false}},
        scales:{x:{title:{display:true,text:xl,font:{size:10},color:'#64748b'},
          ticks:{maxTicksLimit:7,font:{size:8},color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}},
                y:{title:{display:true,text:yl,font:{size:10},color:'#64748b'},
          ticks:{font:{size:8},color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}}}}});
  }
  const Ts=[],hT=[]; for(let t=20;t<=90;t++)      {Ts.push(t);     hT.push(calc(j0,t,P0).h2_kg_hr);}
  const Ps=[],hP=[]; for(let p=1;p<=10;p+=0.1)    {Ps.push(+p.toFixed(1)); hP.push(calc(j0,T0,p).h2_kg_hr);}
  const Js=[],eJ=[]; for(let jj=0.05;jj<=1;jj+=0.01){Js.push(+jj.toFixed(2));eJ.push(calc(jj,T0,P0).eff);}
  mkSC('ch_sT',Ts,hT,'T (°C)','H₂ (kg/hr)','#f59e0b');
  mkSC('ch_sP',Ps,hP,'P (bar)','H₂ (kg/hr)','#8b5cf6');
  mkSC('ch_sj',Js,eJ,'j (A/cm²)','Efficiency (%)','#3b82f6');
}

function phyForecast(r, t_hr){
  const m0=r.h2_kg_hr;
  const lambda=2.5e-5;
  const eps=0.003;
  const omega=2*Math.PI/12;
  const delta=1.2e-4;
  const t_cycle=8;
  const f_T = 1 + 0.0015*(T-60);
  const f_j = 1 - 0.001*Math.max(j-0.2,0);
  return m0 * Math.exp(-lambda*t_hr) * f_T * f_j
       * (1 + eps*Math.sin(omega*t_hr))
       * (1 - delta*Math.floor(t_hr/t_cycle));
}

function updateForecast(r){
  const h0=r.h2_kg_hr;
  const f1=phyForecast(r,1), f6=phyForecast(r,6), f24=phyForecast(r,24);
  st('fc_1h',  f1.toFixed(5)+' kg/hr');
  st('fc_6h',  f6.toFixed(5)+' kg/hr');
  st('fc_24h', f24.toFixed(5)+' kg/hr');
  const d1=(f1-h0)/h0*100, d6=(f6-h0)/h0*100, d24=(f24-h0)/h0*100;
  sh('fc_1t',  `<span style="color:${d1>=0?'#34d399':'#f87171'}">${d1>=0?'+':''}${d1.toFixed(3)}%</span> vs now`);
  sh('fc_6t',  `<span style="color:${d6>=0?'#34d399':'#f87171'}">${d6>=0?'+':''}${d6.toFixed(3)}%</span> vs now`);
  sh('fc_24t', `<span style="color:${d24>=0?'#34d399':'#f87171'}">${d24>=0?'+':''}${d24.toFixed(3)}%</span> vs now`);
  st('fc_params','λ=2.5×10⁻⁵ hr⁻¹ · ε=0.003 · ω=2π/12 rad/hr · δ=1.2×10⁻⁴/cycle');

  const hrs=[],fd=[],fd_hi=[],fd_lo=[];
  for(let h=0;h<=24;h+=0.25){
    hrs.push(h.toFixed(2));
    const fv=phyForecast(r,h), unc=h*0.0003*h0;
    fd.push(fv); fd_hi.push(fv+unc); fd_lo.push(fv-unc);
  }
  const ctx=document.getElementById('ch_fc'); if(!ctx)return;
  if(fcCh) fcCh.destroy();
  fcCh=new Chart(ctx.getContext('2d'),{type:'line',
    data:{labels:hrs, datasets:[
      {label:'Upper',data:fd_hi,borderColor:'rgba(59,130,246,0.2)',backgroundColor:'rgba(59,130,246,0.06)',
       borderWidth:1,pointRadius:0,fill:'+1',tension:0.3,borderDash:[3,3]},
      {label:'H₂ forecast',data:fd,borderColor:'#3b82f6',backgroundColor:'transparent',
       borderWidth:2,pointRadius:0,fill:false,tension:0.3},
      {label:'Lower',data:fd_lo,borderColor:'rgba(59,130,246,0.2)',backgroundColor:'rgba(59,130,246,0.06)',
       borderWidth:1,pointRadius:0,fill:'-1',tension:0.3,borderDash:[3,3]}
    ]},
    options:{responsive:true,maintainAspectRatio:false,
      plugins:{legend:{display:true,labels:{color:'#94a3b8',font:{size:9},boxWidth:10}}},
      scales:{x:{title:{display:true,text:'Time (hr)',font:{size:10},color:'#64748b'},
        ticks:{maxTicksLimit:13,font:{size:8},color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}},
               y:{title:{display:true,text:'H₂ production (kg/hr)',font:{size:10},color:'#64748b'},
        ticks:{font:{size:8},color:'#64748b'},grid:{color:'rgba(255,255,255,0.03)'}}}}});
}

updateUI();
setInterval(()=>{ if(!isPaused) updateUI(); }, 2000);