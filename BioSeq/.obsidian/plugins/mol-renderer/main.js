"use strict";

const { Plugin, Notice } = require("obsidian");

const SVG_NS = "http://www.w3.org/2000/svg";
const SCALE = 58;
const PADDING = 30;
const DEFAULT_BOND_WIDTH = 2.4;
const ATOM_COLORS = {
  C: "#666666",
  H: "#666666",
  "R#": "#666666",
  N: "#3366ee",
  O: "#ff1111",
  F: "#4f8f36",
  Cl: "#13951f",
  Br: "#ee3333",
  I: "#d01ee6",
  S: "#6f7011",
  P: "#c76300",
};

module.exports = class MolRendererPlugin extends Plugin {
  async onload() {
    this.registerMarkdownCodeBlockProcessor("mol", (source, el) => {
      renderMolBlock(source, el);
    });

    this.registerMarkdownCodeBlockProcessor("molfile", (source, el) => {
      renderMolBlock(source, el);
    });

    this.addCommand({
      id: "insert-mol-code-block",
      name: "Insert MOL code block",
      editorCallback: (editor) => {
        editor.replaceSelection("```mol\n\n```\n");
        new Notice("MOL code block inserted");
      },
    });
  }
};

function renderMolBlock(source, containerEl) {
  containerEl.empty();
  containerEl.addClass("mol-renderer");

  try {
    const molecule = parseMolfile(source);
    const svg = buildSvg(molecule);
    const figure = containerEl.createDiv({ cls: "mol-renderer__figure" });
    figure.appendChild(svg);
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    containerEl.createDiv({
      cls: "mol-renderer__error",
      text: `MOL render error: ${message}`,
    });
  }
}

function parseMolfile(source) {
  const lines = source
    .replace(/\r\n?/g, "\n")
    .split("\n")
    .map((line) => line.replace(/\s+$/g, ""));

  while (lines.length > 0 && lines[0].trim() === "") lines.shift();
  while (lines.length > 0 && lines[lines.length - 1].trim() === "") lines.pop();

  const countsIndex = lines.findIndex((line) => parseCountsLine(line) !== null);
  if (countsIndex === -1) {
    throw new Error("could not find a V2000 counts line");
  }

  const counts = parseCountsLine(lines[countsIndex]);
  if (!counts) {
    throw new Error("could not parse atom and bond counts");
  }

  const atomStart = countsIndex + 1;
  const bondStart = atomStart + counts.atomCount;
  if (lines.length < bondStart + counts.bondCount) {
    throw new Error(
      `file ended early: expected ${counts.atomCount} atoms and ${counts.bondCount} bonds`
    );
  }

  const atoms = [];
  for (let index = 0; index < counts.atomCount; index += 1) {
    atoms.push(parseAtomLine(lines[atomStart + index], index + 1));
  }

  const bonds = [];
  for (let index = 0; index < counts.bondCount; index += 1) {
    bonds.push(parseBondLine(lines[bondStart + index], index + 1, atoms.length));
  }

  parseMolProperties(lines.slice(bondStart + counts.bondCount), atoms);

  return { atoms, bonds };
}

function parseCountsLine(line) {
  if (/\bV3000\b/i.test(line)) {
    throw new Error("V3000 molfile is not supported yet");
  }

  const fixedAtomCount = Number.parseInt(line.slice(0, 3), 10);
  const fixedBondCount = Number.parseInt(line.slice(3, 6), 10);
  if (
    Number.isInteger(fixedAtomCount) &&
    Number.isInteger(fixedBondCount) &&
    /\bV2000\b/i.test(line)
  ) {
    return { atomCount: fixedAtomCount, bondCount: fixedBondCount };
  }

  const match = line.match(/^\s*(\d+)\s+(\d+)\b.*\bV2000\b/i);
  if (!match) return null;

  return {
    atomCount: Number.parseInt(match[1], 10),
    bondCount: Number.parseInt(match[2], 10),
  };
}

function parseAtomLine(line, atomNumber) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 4) {
    throw new Error(`atom line ${atomNumber} is too short`);
  }

  const x = Number.parseFloat(parts[0]);
  const y = Number.parseFloat(parts[1]);
  const z = Number.parseFloat(parts[2]);
  const element = parts[3];

  if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
    throw new Error(`atom line ${atomNumber} has invalid coordinates`);
  }

  if (!isSupportedAtomSymbol(element)) {
    throw new Error(`atom line ${atomNumber} has invalid element "${element}"`);
  }

  return { x, y, z, element };
}

function isSupportedAtomSymbol(element) {
  return /^[A-Z][a-z]?$/.test(element) || element === "R#";
}

function parseBondLine(line, bondNumber, atomCount) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 3) {
    throw new Error(`bond line ${bondNumber} is too short`);
  }

  const from = Number.parseInt(parts[0], 10) - 1;
  const to = Number.parseInt(parts[1], 10) - 1;
  const order = Number.parseInt(parts[2], 10);

  if (from < 0 || from >= atomCount || to < 0 || to >= atomCount) {
    throw new Error(`bond line ${bondNumber} references a missing atom`);
  }

  if (![1, 2, 3].includes(order)) {
    throw new Error(`bond line ${bondNumber} has unsupported order ${order}`);
  }

  return { from, to, order };
}

function parseMolProperties(lines, atoms) {
  for (const line of lines) {
    const parts = line.trim().split(/\s+/);
    if (parts.length === 0) continue;
    if (parts[0] !== "M") continue;
    if (parts[1] === "END") break;
    if (parts[1] !== "RGP") continue;

    const pairCount = Number.parseInt(parts[2], 10);
    if (!Number.isInteger(pairCount)) {
      throw new Error("M RGP line has invalid pair count");
    }

    for (let index = 0; index < pairCount; index += 1) {
      const atomNumber = Number.parseInt(parts[3 + index * 2], 10);
      const rGroupNumber = Number.parseInt(parts[4 + index * 2], 10);

      if (!Number.isInteger(atomNumber) || !Number.isInteger(rGroupNumber)) {
        throw new Error("M RGP line has an incomplete atom/group pair");
      }

      const atom = atoms[atomNumber - 1];
      if (!atom) {
        throw new Error(`M RGP line references missing atom ${atomNumber}`);
      }

      atom.rGroupNumber = rGroupNumber;
    }
  }
}

function buildSvg(molecule) {
  if (molecule.atoms.length === 0) {
    throw new Error("molecule has no atoms");
  }

  const points = molecule.atoms.map((atom) => ({
    x: atom.x * SCALE,
    y: -atom.y * SCALE,
  }));

  const bounds = getBounds(points);
  const width = Math.max(120, bounds.maxX - bounds.minX + PADDING * 2);
  const height = Math.max(120, bounds.maxY - bounds.minY + PADDING * 2);
  const normalized = points.map((point) => ({
    x: point.x - bounds.minX + PADDING,
    y: point.y - bounds.minY + PADDING,
  }));

  const svg = createSvgElement("svg");
  svg.setAttribute("xmlns", SVG_NS);
  svg.setAttribute("viewBox", `0 0 ${round(width)} ${round(height)}`);
  svg.setAttribute("width", round(width));
  svg.setAttribute("height", round(height));
  svg.setAttribute("style", `width: min(100%, ${round(width)}px); height: auto; color: var(--text-normal, #222);`);
  svg.setAttribute("role", "img");
  svg.setAttribute("aria-label", `Molecular structure with ${molecule.atoms.length} atoms and ${molecule.bonds.length} bonds`);

  const bondsGroup = createSvgElement("g");
  const labelsGroup = createSvgElement("g");
  svg.appendChild(bondsGroup);
  svg.appendChild(labelsGroup);

  for (const bond of molecule.bonds) {
    drawBond(bondsGroup, normalized[bond.from], normalized[bond.to], bond.order);
  }

  molecule.atoms.forEach((atom, index) => {
    drawAtomLabel(labelsGroup, normalized[index], atom, index, molecule.bonds);
  });

  return svg;
}

function getBounds(points) {
  return points.reduce(
    (bounds, point) => ({
      minX: Math.min(bounds.minX, point.x),
      maxX: Math.max(bounds.maxX, point.x),
      minY: Math.min(bounds.minY, point.y),
      maxY: Math.max(bounds.maxY, point.y),
    }),
    {
      minX: Number.POSITIVE_INFINITY,
      maxX: Number.NEGATIVE_INFINITY,
      minY: Number.POSITIVE_INFINITY,
      maxY: Number.NEGATIVE_INFINITY,
    }
  );
}

function drawBond(group, from, to, order) {
  if (order === 1) {
    appendBondLine(group, from, to, 0);
    return;
  }

  if (order === 2) {
    appendBondLine(group, from, to, -4);
    appendBondLine(group, from, to, 4);
    return;
  }

  appendBondLine(group, from, to, -6);
  appendBondLine(group, from, to, 0);
  appendBondLine(group, from, to, 6);
}

function appendBondLine(group, from, to, offset) {
  const shifted = offsetLine(from, to, offset);
  const line = createSvgElement("line");
  line.classList.add("mol-renderer__bond");
  line.setAttribute("x1", round(shifted.from.x));
  line.setAttribute("y1", round(shifted.from.y));
  line.setAttribute("x2", round(shifted.to.x));
  line.setAttribute("y2", round(shifted.to.y));
  line.setAttribute("stroke", "currentColor");
  line.setAttribute("stroke-width", String(DEFAULT_BOND_WIDTH));
  line.setAttribute("stroke-linecap", "round");
  group.appendChild(line);
}

function offsetLine(from, to, offset) {
  const dx = to.x - from.x;
  const dy = to.y - from.y;
  const length = Math.hypot(dx, dy) || 1;
  const nx = (-dy / length) * offset;
  const ny = (dx / length) * offset;

  return {
    from: { x: from.x + nx, y: from.y + ny },
    to: { x: to.x + nx, y: to.y + ny },
  };
}

function drawAtomLabel(group, point, atom, atomIndex, bonds) {
  if (!shouldDrawLabel(atom, atomIndex, bonds)) return;
  const labelParts = getAtomLabelParts(atom, atomIndex, bonds);

  const background = createSvgElement("text");
  background.classList.add("mol-renderer__atom-bg");
  background.setAttribute("x", round(point.x));
  background.setAttribute("y", round(point.y));
  background.setAttribute("fill", "var(--background-primary, #fff)");
  background.setAttribute("stroke", "var(--background-primary, #fff)");
  background.setAttribute("stroke-width", "6");
  background.setAttribute("stroke-linejoin", "round");
  background.setAttribute("font-family", "Arial, sans-serif");
  background.setAttribute("font-size", "16");
  background.setAttribute("font-weight", "600");
  background.setAttribute("dominant-baseline", "central");
  background.setAttribute("text-anchor", "middle");
  appendLabelParts(background, labelParts);
  group.appendChild(background);

  const label = createSvgElement("text");
  label.classList.add("mol-renderer__atom-label");
  label.classList.add(`mol-renderer__atom-label--${atom.element}`);
  label.setAttribute("x", round(point.x));
  label.setAttribute("y", round(point.y));
  const color = ATOM_COLORS[atom.element] ?? "currentColor";
  label.setAttribute("fill", color);
  label.setAttribute("style", `fill: ${color};`);
  label.setAttribute("font-family", "Arial, sans-serif");
  label.setAttribute("font-size", "16");
  label.setAttribute("font-weight", "600");
  label.setAttribute("dominant-baseline", "central");
  label.setAttribute("text-anchor", "middle");
  appendLabelParts(label, labelParts);
  group.appendChild(label);
}

function getAtomLabelParts(atom, atomIndex, bonds) {
  if (atom.element === "R#" && Number.isInteger(atom.rGroupNumber)) {
    return [{ text: `R${atom.rGroupNumber}` }];
  }

  if (atom.element === "R#") {
    return [{ text: "R" }];
  }

  if (atom.element === "C") {
    const valenceUsed = getBondOrderSum(atomIndex, bonds);
    const hydrogenCount = Math.max(0, 4 - valenceUsed);

    if (hydrogenCount === 0) {
      return [{ text: "C" }];
    }

    return [
      { text: "C" },
      { text: "H" },
      ...(hydrogenCount > 1 ? [{ text: String(hydrogenCount), subscript: true }] : []),
    ];
  }

  return [{ text: atom.element }];
}

function getAtomLabel(atom) {
  return getAtomLabelParts(atom, -1, [])
    .map((part) => part.text)
    .join("");
}

function appendLabelParts(textElement, parts) {
  textElement.textContent = "";

  for (const part of parts) {
    const tspan = createSvgElement("tspan");
    tspan.textContent = part.text;

    if (part.subscript) {
      tspan.setAttribute("font-size", "11");
      tspan.setAttribute("baseline-shift", "sub");
    }

    textElement.appendChild(tspan);
  }
}

function getBondOrderSum(atomIndex, bonds) {
  return bonds.reduce((sum, bond) => {
    if (bond.from !== atomIndex && bond.to !== atomIndex) {
      return sum;
    }

    return sum + bond.order;
  }, 0);
}

function shouldDrawLabel(atom, atomIndex, bonds) {
  return true;
}

function createSvgElement(tagName) {
  return document.createElementNS(SVG_NS, tagName);
}

function round(value) {
  return String(Math.round(value * 100) / 100);
}
