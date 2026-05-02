/* FRNC-5Pro  –  main.js
   Handles dynamic convection bank UI and demo data loading.
*/

let bankCount = 0;

function addBank() {
  bankCount++;
  const idx = bankCount - 1;
  const container = document.getElementById('conv-banks-container');
  const div = document.createElement('div');
  div.className = 'conv-bank-card';
  div.id = `bank-card-${idx}`;

  div.innerHTML = `
    <div class="conv-bank-header">
      管排 ${bankCount} / Bank ${bankCount}
      <button type="button" class="remove-bank-btn" onclick="removeBank(${idx})">
        ✕ 删除 / Remove
      </button>
    </div>
    <div class="grid-3">
      <div class="field">
        <label>外径 OD <span class="unit">(mm)</span></label>
        <input type="number" name="conv_bank_${idx}_tube_od_mm" value="114.3" step="0.1" min="10"/>
      </div>
      <div class="field">
        <label>壁厚 WT <span class="unit">(mm)</span></label>
        <input type="number" name="conv_bank_${idx}_tube_wt_mm" value="8.0" step="0.5"/>
      </div>
      <div class="field">
        <label>横向节距 S_T <span class="unit">(mm)</span></label>
        <input type="number" name="conv_bank_${idx}_pitch_trans_mm" value="220" step="5"/>
      </div>
      <div class="field">
        <label>纵向节距 S_L <span class="unit">(mm)</span></label>
        <input type="number" name="conv_bank_${idx}_pitch_long_mm" value="220" step="5"/>
      </div>
      <div class="field">
        <label>有效管长 <span class="unit">(m)</span></label>
        <input type="number" name="conv_bank_${idx}_tube_length_m" value="10" step="0.5"/>
      </div>
      <div class="field">
        <label>每排管数 / Tubes per Row</label>
        <input type="number" name="conv_bank_${idx}_n_tubes_per_row" value="12" step="1" min="1"/>
      </div>
      <div class="field">
        <label>排数 / No. Rows</label>
        <input type="number" name="conv_bank_${idx}_n_rows" value="4" step="1" min="1"/>
      </div>
      <div class="field">
        <label>排列方式 / Arrangement</label>
        <select name="conv_bank_${idx}_arrangement">
          <option value="staggered" selected>交叉排列 / Staggered</option>
          <option value="inline">顺排 / Inline</option>
        </select>
      </div>
      <div class="field">
        <label>管材 / Tube Material</label>
        <select name="conv_bank_${idx}_tube_material">
          <option value="carbon_steel" selected>碳钢 / Carbon Steel</option>
          <option value="Cr5Mo">5Cr½Mo (P5)</option>
          <option value="SS304">304 Stainless</option>
          <option value="SS316">316 Stainless</option>
        </select>
      </div>
    </div>
    <details style="margin-top:10px">
      <summary style="cursor:pointer;font-size:12px;color:#666">
        ▶ 翅片参数 / Fin Parameters (可选 / optional)
      </summary>
      <div class="grid-3" style="margin-top:8px">
        <div class="field">
          <label>翅片高度 / Fin Height <span class="unit">(mm, 0=裸管)</span></label>
          <input type="number" name="conv_bank_${idx}_fin_height_mm" value="0" step="1" min="0"/>
        </div>
        <div class="field">
          <label>翅片厚度 / Fin Thickness <span class="unit">(mm)</span></label>
          <input type="number" name="conv_bank_${idx}_fin_thickness_mm" value="2.5" step="0.5"/>
        </div>
        <div class="field">
          <label>翅片密度 / Fin Density <span class="unit">(fins/m)</span></label>
          <input type="number" name="conv_bank_${idx}_fin_pitch_per_m" value="197" step="10"/>
        </div>
      </div>
    </details>
  `;

  container.appendChild(div);
}


function removeBank(idx) {
  const el = document.getElementById(`bank-card-${idx}`);
  if (el) el.remove();
}


function loadDemo() {
  /* Fill in a typical crude preheat furnace demo case */
  const set = (name, val) => {
    const el = document.querySelector(`[name="${name}"]`);
    if (el) el.value = val;
  };

  set('case_name', 'Demo – 常压蒸馏加热炉');
  set('heater_service', '常压炉 / Atmospheric Crude Heater');
  set('fuel_type', 'natural_gas');
  set('excess_air_pct', 15);
  set('W_fuel_kg_s', 0.28);
  set('T_air_C', 15);
  set('loss_fraction_pct', 2);
  set('firebox_type', 'box');
  set('firebox_length_m', 12);
  set('firebox_width_m', 5);
  set('firebox_height_m', 7);
  set('n_tubes_radiant', 40);
  set('rad_tube_od_mm', 168.3);
  set('rad_tube_wt_mm', 9);
  set('rad_tube_pitch_mm', 336.6);
  set('n_radiant_rows', 1);
  set('rad_tube_length_m', 12);
  set('rad_tube_material', 'Cr5Mo');
  set('process_fluid', 'crude_oil');
  set('W_process_kg_s', 80);
  set('T_proc_in_C', 120);
  set('T_proc_out_C', '');
  set('peak_flux_factor', 1.8);

  /* Add one convection bank */
  document.getElementById('conv-banks-container').innerHTML = '';
  bankCount = 0;
  addBank();
  // Set convection bank values
  set('conv_bank_0_tube_od_mm', 114.3);
  set('conv_bank_0_tube_wt_mm', 8);
  set('conv_bank_0_pitch_trans_mm', 220);
  set('conv_bank_0_pitch_long_mm', 220);
  set('conv_bank_0_tube_length_m', 12);
  set('conv_bank_0_n_tubes_per_row', 14);
  set('conv_bank_0_n_rows', 6);

  alert('示例参数已加载！\nDemo parameters loaded.');
}
