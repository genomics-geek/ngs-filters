const { allelic_balance_high_quality } = require('./index')
const { high_quality } = require('./index')
const { uniparental_disomy } = require('./index')
const { denovo } = require('./index')
const { x_linked_denovo } = require('./index')
const { homozygous_recessive } = require('./index')
const { x_linked_homozygous_recessive } = require('./index')
const { compound_heterozygous_side } = require('./index')
const { remove_compound_heterozygous_side } = require('./index')
const { in_hgmd } = require('./index')
const { nonsynonymous } = require('./index')
const { cohort_af } = require('./index')

test.each([
  ['allelic_balance_high_quality: homozygous REF with high AB', { alts: 0, AB: 0.2 }, false],
  ['allelic_balance_high_quality: homozygous REF with appropriate AB', { alts: 0, AB: 0.001 }, true],
  ['allelic_balance_high_quality: heterozygous with low AB', { alts: 1, AB: 0.01 }, false],
  ['allelic_balance_high_quality: heterozygous with high AB', { alts: 1, AB: 0.9 }, false],
  ['allelic_balance_high_quality: heterozygous with appropriate AB', { alts: 1, AB: 0.3 }, true],
  ['allelic_balance_high_quality: homozygous REF with appropriate AB', { alts: 2, AB: 0.99 }, true],
])('%s: %j equals %p', (a, b, expected) => {
  expect(allelic_balance_high_quality(b)).toBe(expected)
})

test.each([
  ['high_quality: high DP but low GQ', { DP: 100, GQ: 5 }, false],
  ['high_quality: low DP but high GQ', { DP: 1, GQ: 75 }, false],
])('%s: %j equals %p', (a, b, expected) => {
  expect(high_quality(b)).toBe(expected)
})

test.each([
  ['uniparental_disomy: can handle missing dad', { alts: 2 }, { alts: 1 }, undefined, false],
  ['uniparental_disomy: can handle missing mom', { alts: 2 }, undefined, { alts: 1 }, false],
  ['uniparental_disomy: Not inherited (homozygous recessive)', { alts: 2 }, { alts: 1 }, { alts: 1 }, false],
  ['uniparental_disomy: 2 ALT alleles inherited from mom', { alts: 2 }, { alts: 2 }, { alts: 0 }, true],
  ['uniparental_disomy: 2 ALT alleles inherited from dad', { alts: 2 }, { alts: 0 }, { alts: 2 }, true],
  ['uniparental_disomy: 2 REF alleles inherited from mom', { alts: 0 }, { alts: 0 }, { alts: 2 }, true],
  ['uniparental_disomy: 2 REF alleles inherited from dad', { alts: 0 }, { alts: 2 }, { alts: 0 }, true],
  ['uniparental_disomy: proband is heterozygous', { alts: 1 }, { alts: 2 }, { alts: 0 }, false],
  [
    'uniparental_disomy: proband is reference and parents are heterozygous',
    { alts: 0 },
    { alts: 1 },
    { alts: 1 },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(uniparental_disomy(kid, mom, dad)).toBe(expected)
})

test.each([
  ['denovo: can handle missing parents', { alts: 1 }, { alts: 0 }, undefined, true],
  ['denovo: can handle missing parents', { alts: 1 }, undefined, undefined, true],
  ['denovo: parents are both homozygous REF', { alts: 1 }, { alts: 0 }, { alts: 0 }, true],
  ['denovo: proband is homozygous', { alts: 2 }, { alts: 0 }, { alts: 0 }, false],
  ['denovo: mom is heterozygous', { alts: 1 }, { alts: 1 }, { alts: 0 }, false],
  ['denovo: dad is heterozygous', { alts: 1 }, { alts: 0 }, { alts: 1 }, false],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(denovo(kid, mom, dad)).toBe(expected)
})

test.each([
  ['x_linked_denovo: can handle missing parents', { alts: 1 }, { alts: 0 }, undefined, true],
  ['x_linked_denovo: can handle missing parents', { alts: 1 }, undefined, undefined, true],
  ['x_linked_denovo: can handle missing sex', { alts: 1 }, { alts: 0 }, { alts: 0 }, true],
  ['x_linked_denovo: parents are both homozygous REF', { alts: 1 }, { alts: 0 }, { alts: 0 }, true],
  ['x_linked_denovo: proband is homozygous', { alts: 2 }, { alts: 0 }, { alts: 0 }, true],
  ['x_linked_denovo: proband is female', { alts: 1, sex: 'female' }, { alts: 0 }, { alts: 0 }, false],
  ['x_linked_denovo: mom is heterozygous', { alts: 1, sex: 'male' }, { alts: 1 }, { alts: 0 }, false],
  ['x_linked_denovo: dad is hemizygous', { alts: 1, sex: 'male' }, { alts: 0 }, { alts: 1 }, false],
  ['x_linked_denovo: proband is homozygous REF', { alts: 0, sex: 'male' }, { alts: 0 }, { alts: 1 }, false],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(x_linked_denovo(kid, mom, dad)).toBe(expected)
})

test.each([
  ['homozygous_recessive: can handle missing parents', { alts: 2 }, undefined, undefined, true],
  ['homozygous_recessive: mom is heterozygous and dad is heterozygous', { alts: 2 }, { alts: 1 }, { alts: 1 }, true],
  ['homozygous_recessive: mom is homozygous and dad is heterozygous', { alts: 2 }, { alts: 2 }, { alts: 1 }, true],
  ['homozygous_recessive: mom is heterozygous and dad is homozygous', { alts: 2 }, { alts: 1 }, { alts: 2 }, true],
  ['homozygous_recessive: proband is heterozygous', { alts: 1 }, { alts: 1 }, { alts: 2 }, false],
  ['homozygous_recessive: proband is reference', { alts: 0 }, { alts: 1 }, { alts: 2 }, false],
  [
    'homozygous_recessive: proband is homozygous but mom is homozygous REF',
    { alts: 2 },
    { alts: 0 },
    { alts: 2 },
    false,
  ],
  [
    'homozygous_recessive: proband is homozygous but dad is homozygous REF',
    { alts: 2 },
    { alts: 1 },
    { alts: 0 },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(homozygous_recessive(kid, mom, dad)).toBe(expected)
})

test.each([
  ['x_linked_homozygous_recessive: can handle missing parents', { alts: 1 }, { alts: 1 }, undefined, true],
  ['x_linked_homozygous_recessive: can handle missing parents', { alts: 1 }, undefined, undefined, true],
  ['x_linked_homozygous_recessive: can handle missing sex', { alts: 1 }, { alts: 1 }, { alts: 0 }, true],
  ['x_linked_homozygous_recessive: parents are both homozygous REF', { alts: 1 }, { alts: 0 }, { alts: 0 }, false],
  [
    'x_linked_homozygous_recessive: proband is homozygous and mom is heterozygous',
    { alts: 2 },
    { alts: 1 },
    { alts: 0 },
    true,
  ],
  [
    'x_linked_homozygous_recessive: proband is homozygous, mom is heterozygous, and dad is heterozygous',
    { alts: 2 },
    { alts: 1 },
    { alts: 1 },
    false,
  ],
  ['x_linked_homozygous_recessive: proband is female', { alts: 1, sex: 'female' }, { alts: 0 }, { alts: 0 }, false],
  [
    'x_linked_homozygous_recessive: proband is homozygous REF',
    { alts: 0, sex: 'male' },
    { alts: 0 },
    { alts: 0 },
    false,
  ],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(x_linked_homozygous_recessive(kid, mom, dad)).toBe(expected)
})

test.each([
  ['compound_heterozygous_side: can handle missing parents', { alts: 1 }, undefined, undefined, true],
  [
    'compound_heterozygous_side: mom is heterozygous and dad is heterozygous',
    { alts: 1 },
    { alts: 1 },
    { alts: 1 },
    true,
  ],
  [
    'compound_heterozygous_side: mom is homozygous and dad is heterozygous',
    { alts: 1 },
    { alts: 2 },
    { alts: 1 },
    true,
  ],
  [
    'compound_heterozygous_side: mom is heterozygous and dad is homozygous',
    { alts: 1 },
    { alts: 1 },
    { alts: 2 },
    true,
  ],
  ['compound_heterozygous_side: mom and dad are both homozygous REF', { alts: 1 }, { alts: 0 }, { alts: 0 }, false],
  [
    'compound_heterozygous_side: mom is undefined and dad is homozygous REF',
    { alts: 1 },
    undefined,
    { alts: 0 },
    false,
  ],
  [
    'compound_heterozygous_side: dad is undefined and mom is homozygous REF',
    { alts: 1 },
    { alts: 0 },
    undefined,
    false,
  ],
  ['compound_heterozygous_side: proband is not heterozygous', { alts: 0 }, { alts: 0 }, undefined, false],
])('%s: Kid: %j Mom: %j Dad: %j equals %p', (a, kid, mom, dad, expected) => {
  expect(compound_heterozygous_side(kid, mom, dad)).toBe(expected)
})

test.each([
  ['remove_compound_heterozygous_side: can handle defaults', { compound_heterozygous_side: 'id' }, false],
  [
    'remove_compound_heterozygous_side: will remove compound_heterozygous_side',
    { compound_heterozygous_side: 'id' },
    false,
  ],
  ['remove_compound_heterozygous_side: will not remove compound_heterozygous_side', {}, true],
])('%s: INFO: %j equals %p', (a, INFO, expected) => {
  expect(remove_compound_heterozygous_side(INFO)).toBe(expected)
})

test.each([
  ['in_hgmd: can handle missing CSQ', {}, false],
  ['in_hgmd: HGMD ID is not set', { CSQ: 'x||' }, false],
  ['in_hgmd: HGMD ID is set', { CSQ: 'x|x|' }, true],
])('%s: %j equals %p', (a, info, expected) => {
  expect(in_hgmd(info, 'CSQ', 1)).toBe(expected)
})

test.each([
  ['in_hgmd: can handle default parameters', { CSQ: '||||||||||||||||||||||||||||x||||||||' }, true],
  ['in_hgmd: can handle default parameters', { CSQ: '||||||||||||||||||||||||||||||||||||' }, false],
  ['in_hgmd: can handle default parameters', { CSQ: '||||' }, false],
])('%s: %j equals %p', (a, info, expected) => {
  expect(in_hgmd(info)).toBe(expected)
})

test.each([
  ['nonsynonymous: can handle missing CSQ', {}, false],
  ['nonsynonymous: can filter out variants that are not nonsynonymous', { CSQ: 'x|synonymous|' }, false],
  ['nonsynonymous: will retain nonsynonymous filters', { CSQ: 'x|splice|' }, true],
])('%s: %j equals %p', (a, info, expected) => {
  expect(nonsynonymous(info)).toBe(expected)
})

test.each([
  ['cohort_af: can handle missing CSQ', {}, false],
  ['cohort_af: AF is not set', { CSQ: 'x||' }, true],
  ['cohort_af: AF is below cutoff', { CSQ: 'x|0.001|' }, true],
  ['cohort_af: AF is above cutoff', { CSQ: 'x|0.5|' }, false],
])('%s: %j equals %p', (a, info, expected) => {
  expect(cohort_af(info, 'CSQ', 1)).toBe(expected)
})

test.each([
  ['cohort_af: can handle default parameters', { CSQ: '|||||||||||||||||||||||||||||||||||0.2|' }, false],
  ['cohort_af: can handle default parameters', { CSQ: '|||||||||||||||||||||||||||||||||||0.005|' }, true],
  ['cohort_af: can handle default parameters', { CSQ: '||||||||||||||||||||||||||||||||||||' }, true],
  ['cohort_af: can handle default parameters', { CSQ: '||||||||||||' }, true],
])('%s: %j equals %p', (a, info, expected) => {
  expect(cohort_af(info)).toBe(expected)
})
