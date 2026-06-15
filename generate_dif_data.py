import numpy as np
import pandas as pd


def simulate_rasch_dif(
    n_students=10000,
    n_items=60,
    dif_items_moderate=6,
    dif_items_large=6,
    # Taken from winsteps
    moderate_range=(0.43, 0.64),
    large_range=(0.64, 0.90),
    seed=67
):
    rng = np.random.default_rng(seed)

    if dif_items_moderate + dif_items_large > n_items:
        raise ValueError("Total DIF items cannot exceed n_items")

    student_id = np.arange(1, n_students + 1)
    group = rng.integers(0, 2, size=n_students)  # 0/1 groups
    theta = rng.normal(0, 1, size=n_students)

    item_ids = [f"Item_{j+1:02d}" for j in range(n_items)]
    b = rng.normal(0, 1, size=n_items)

    n_dif = dif_items_moderate + dif_items_large
    dif_idx = rng.choice(n_items, size=n_dif, replace=False)

    moderate_idx = dif_idx[:dif_items_moderate]
    large_idx = dif_idx[dif_items_moderate:]

    delta = np.zeros(n_items)

    if dif_items_moderate > 0:
        delta[moderate_idx] = rng.uniform(
            moderate_range[0], moderate_range[1], size=dif_items_moderate
        )

    if dif_items_large > 0:
        delta[large_idx] = rng.uniform(
            large_range[0], large_range[1], size=dif_items_large
        )

    # Group 1 sees DIF items as harder by delta_j
    # Group 0 sees baseline difficulty b_j
    responses = np.zeros((n_students, n_items), dtype=int)

    for j in range(n_items):
        b_j_group = np.where(group == 1, b[j] + delta[j], b[j])
        p = 1.0 / (1.0 + np.exp(-(theta - b_j_group)))
        responses[:, j] = rng.binomial(1, p, size=n_students)

    response_df = pd.DataFrame(responses, columns=item_ids)
    response_df.insert(0, "group", group)
    response_df.insert(0, "student_id", student_id)

    truth_df = pd.DataFrame({
        "item_id": item_ids,
        "b_base": b,
        "delta_group1": delta,
        "dif_flag": delta > 0,
        "dif_type": np.where(
            delta >= large_range[0], "***",
            np.where(delta >= moderate_range[0], "**", "")
        )
    })

    latent_df = pd.DataFrame({
        "student_id": student_id,
        "group": group,
        "theta": theta
    })

    return response_df, truth_df, latent_df


if __name__ == "__main__":
    response_df, truth_df, latent_df = simulate_rasch_dif(
        n_students=10000,
        n_items=60,
        dif_items_moderate=6,
        dif_items_large=6,
        seed=67
    )

    response_df.to_csv("simulated_rasch_dif_responses.csv", index=False)
    truth_df.to_csv("simulated_rasch_dif_truth.csv", index=False)
    latent_df.to_csv("simulated_rasch_dif_latent.csv", index=False)

    print(response_df.head())
    print()
    print(truth_df)