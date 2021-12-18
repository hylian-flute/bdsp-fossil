const WALL_WIDTH: usize = 13;
const WALL_HEIGHT: usize = 10;
const WALL_HIT_POINT: u64 = 50;
const GENE_LENGTH: usize = 50;

#[derive(PartialEq, Clone)]
enum Tool {
  Pickel,
  Hammer,
}

#[derive(Clone)]
struct Action {
  tool: Tool,
  point: usize,
}

impl Action {
  pub fn new() -> Action {
    Action {
      tool: if rand::random::<bool>() { Tool::Pickel } else { Tool::Hammer },
      point: ((WALL_WIDTH * WALL_HEIGHT) as f64 * rand::random::<f64>()) as usize,
    }
  }
}

struct Wall {
  blocks: [u64; WALL_WIDTH * WALL_HEIGHT],
  hit_point: u64,
}

impl Wall {
  fn shave(&mut self, point: usize, x: isize, y: isize, depth: u64) {
    let point = point as isize;
    let wall_width = WALL_WIDTH as isize;
    if point + x < 0 || point / wall_width != (point + x) / wall_width {
      return;
    }
    let target_point = point + x + wall_width * y;
    if target_point < 0 || target_point >= wall_width * WALL_HEIGHT as isize {
      return;
    }
    let target_point = target_point as usize;
    self.blocks[target_point] = match self.blocks[target_point].checked_sub(depth) {
      Some(v) => v,
      None => 0,
    };
  }
  pub fn act(&mut self, action: &Action) {
    if self.hit_point <= 0 { return; };
    let damage: u64; 
    match action.tool {
      Tool::Pickel => {
        damage = 1;
        self.shave(action.point, 0, -1, 1);
        self.shave(action.point, -1, 0, 1);
        self.shave(action.point, 0, 0, 2);
        self.shave(action.point, 1, 0, 1);
        self.shave(action.point, 0, 1, 1);
      },
      Tool::Hammer => {
        damage = 2;
        self.shave(action.point, -1, -1, 1);
        self.shave(action.point, 0, -1, 2);
        self.shave(action.point, 1, -1, 1);
        self.shave(action.point, -1, 0, 2);
        self.shave(action.point, 0, 0, 2);
        self.shave(action.point, 1, 0, 2);
        self.shave(action.point, -1, 1, 1);
        self.shave(action.point, 0, 1, 2);
        self.shave(action.point, 1, 1, 1);
      },
    };
    self.hit_point = match self.hit_point.checked_sub(damage) {
      Some(v) => v,
      None => 0,
    };
  }
  pub fn evaluate(&self) -> u64 {
    let mut result = 0u64;
    for x in 0..WALL_WIDTH - 1 {
      for y in 0..WALL_HEIGHT - 1 {
        let in_range_blocks = [
          self.blocks[WALL_WIDTH * y + x],
          self.blocks[WALL_WIDTH * y + x + 1],
          self.blocks[WALL_WIDTH * (y + 1) + x],
          self.blocks[WALL_WIDTH * (y + 1) + x + 1],
        ];
        let min = in_range_blocks.iter().fold(u64::MAX, |min, v| if *v < min { *v } else { min });
        if min < 2 { result += 1; };
      }
    }
    result
  }
}

#[derive(Clone)]
struct Gene {
  actions: Vec<Action>,
}

impl Gene {
  pub fn new() -> Gene {
    let mut actions = Vec::with_capacity(GENE_LENGTH);
    for _ in 0..GENE_LENGTH { actions.push(Action::new()) }
    Gene {
      actions: actions,
    }
  }
  pub fn evaluate(&self) -> u64 {
    let mut wall = Wall {
      blocks: [
        6, 6, 4, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6,
        6, 4, 4, 4, 2, 2, 2, 2, 4, 4, 6, 6, 6,
        4, 4, 4, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4,
        2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4,
        2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4,
        2, 4, 4, 4, 4, 2, 2, 2, 2, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 4, 4, 4,
        6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 2, 2, 2,
        6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 2, 2,
        6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 2, 2, 2,
      ],
      hit_point: WALL_HIT_POINT,
    };
    for action in self.actions.iter() {
      wall.act(action);
    }
    wall.evaluate()
  }
  pub fn cross(a: &Gene, b: &Gene) -> Gene {
    if a.actions.len() != b.actions.len() {
      panic!("Crossing of genes of different lengths.")
    }
    let mut actions = Vec::with_capacity(a.actions.len());
    for i in 0..a.actions.len() {
      actions.push(if rand::random::<f64>() < 0.5 {
        a.actions[i].clone()
      } else {
        b.actions[i].clone()
      });
    }
    Gene {
      actions: actions,
    }
  }
  pub fn mutation(&mut self) {
    let index = (self.actions.len() as f64 * rand::random::<f64>()) as usize;
    self.actions[index] = Action::new();
  }
}

struct GeneticAlgorithm {
  genes_num: usize,
  generations_num: usize,
  mutation_probability: f64,
}

impl GeneticAlgorithm {
  pub fn run(&self) -> (Gene, u64) {
    let mut genes: Vec<Gene> = Vec::with_capacity(self.genes_num);
    for _ in 0..self.genes_num { genes.push(Gene::new()) }
    let mut best_gene = genes[0].clone();
    let mut best_gene_value = 0u64;
    for _ in 0..self.generations_num {
      let values: Vec<u64> = genes.iter().map(|gene| gene.evaluate()).collect();
      let (min, sum) = values.iter().fold(
        (u64::MAX, 0),
        |acc, v| (if *v < acc.0 { *v } else { acc.0 }, acc.1 + *v)
      );

      let sum_of_diffs_from_min = sum - self.genes_num as u64 * min;
      let possibilities: Vec<f64> = match sum_of_diffs_from_min {
        0 => vec![1f64 / self.genes_num as f64; self.genes_num],
        _ => values
          .iter()
          .map(|value| (value - min) as f64 / sum_of_diffs_from_min as f64)
          .collect(),
      };

      let mut new_genes = Vec::with_capacity(self.genes_num);
      for _ in 0..self.genes_num {
        let mut parents = Vec::with_capacity(2);
        for _ in 0..2 {
          let random = rand::random::<f64>();
          let mut cumulative_possibility = 0f64;
          for i in 0..self.genes_num {
            cumulative_possibility += possibilities[i];
            if cumulative_possibility >= random {
              parents.push(&genes[i]);
              break;
            }
          }
        }
        let mut gene = Gene::cross(parents[0], parents[1]);
        if rand::random::<f64>() < self.mutation_probability {
          gene.mutation();
        }
        new_genes.push(gene);
      }
      genes = new_genes;

      let mut current_best_gene = &genes[0];
      let mut current_best_gene_value = 0u64;
      for i in 0..self.genes_num {
        if values[i] >= current_best_gene_value {
          current_best_gene_value = values[i];
          current_best_gene = &genes[i];
        }
      }
      if current_best_gene_value >= best_gene_value {
        best_gene_value = current_best_gene_value;
        best_gene = current_best_gene.clone();
      }
    }

    (best_gene, best_gene_value)
  }
}

fn main() {
  let genetic_algorithm = GeneticAlgorithm {
    genes_num: 32,
    generations_num: 100000,
    mutation_probability: 1.0 / 128.0,
  };
  let (gene, value) = genetic_algorithm.run();
  println!("gene:");
  for action in gene.actions.iter() {
    println!(
      "  {}, ({}, {})",
      match action.tool {
        Tool::Pickel => "Pickel",
        Tool::Hammer => "Hammer",
      },
      action.point % WALL_WIDTH,
      action.point / WALL_WIDTH
    );
  }
  println!("value: {}", &value);
}
